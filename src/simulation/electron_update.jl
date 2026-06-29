# update useful quantities relevant for potential, electron energy and fluid solve
function update_electrons!(params, config, t = 0)
    (; Tev, ne, nϵ, νan, νc, νen, νei, radial_loss_frequency, Z_eff, νiz, νex, νew_momentum, κ) = params.cache
    (; source_energy, wall_loss_model, conductivity_model, anom_model) = config

    t_ms = t*1e3
    @printf("updating electrons at time: %.5f ms\n", t_ms)

    # Update electron temperature given new density
    update_temperature!(Tev, nϵ, ne)

    # Update collision frequencies
    if (params.electron_ion_collisions)
        freq_electron_ion!(νei, ne, Tev, Z_eff)
    end

    # Add electron-neutral MEX collisions
    νen .= 0
    for (i, coll) in zip(params.electron_neutral_indices, params.electron_neutral_collisions)
        neutral = params.fluid_array[i]
        freq_electron_neutral!(νen, coll, neutral, Tev)
    end

    freq_electron_classical!(νc, νen, νei, νiz, νex, params.landmark)
    freq_electron_wall!(radial_loss_frequency, νew_momentum, wall_loss_model, params)

    # Update anomalous transport
    t > 0 && anom_model(νan, params, config)

    # Update mobility, discharge current, potential, and electron velocity
    update_electrical_vars!(params)

    # update the thermal conductivity
    conductivity_model(κ, params)

    # Update the electron energy density,  temperature and pressure
    update_electron_energy!(params, wall_loss_model, source_energy, params.dt[])

    return
end

function freq_electron_wall!(radial_loss_frequency, νew_momentum, wall_loss_model, params)
    (; thruster, grid, transition_length) = params
    L_ch = thruster.geometry.channel_length
    # Update wall collisions
    return @inbounds for i in eachindex(radial_loss_frequency)
        # Compute wall collision frequency, with transition function to force no momentum wall collisions in plume
        radial_loss_frequency[i] = freq_electron_wall(
            wall_loss_model, params, i,
        )
        νew_momentum[i] = radial_loss_frequency[i] * linear_transition(
            grid.cell_centers[i], L_ch, transition_length, 1.0, 0.0,
        )
    end
end

function update_electrical_vars!(params)
    (; cache, anom_smoothing_iters, landmark, grid) = params
    (;
        cell_cache_1, νan, νe, νc, μ, B, νew_momentum, anom_multiplier,
        Vs, ue, ji, channel_area, ne, Id, Id_L_IE, Id_IE_R, K, pe, ∇pe, ϕ, ∇ϕ,
    ) = cache

    # Smooth anomalous transport model
    if anom_smoothing_iters > 0
        smooth!(νan, cell_cache_1, iters = anom_smoothing_iters)
    end

    # Update the anomalous collision frequency multiplier to match target current
    anom_multiplier[] = exp(
        apply_controller(
            params.simulation.current_control, Id[], log(anom_multiplier[]), params.dt[],
        )
    )

    @inbounds for i in eachindex(νan)
        # Multiply by anom anom multiplier for PID control
        νan[i] *= anom_multiplier[]

        # Compute total collision frequency and electron mobility
        νe[i] = νc[i] + νan[i] + νew_momentum[i]
        μ[i] = electron_mobility(νe[i], B[i])
    end

    # Compute anode sheath potential
    # TODO: should this go here?
    Vs[] = anode_sheath_potential(params)

    # Compute the discharge current by integrating the momentum equation over the whole domain
    V_L = params.discharge_voltage + Vs[]
    V_R = params.cathode_coupling_voltage # usually zero

    apply_drag = !landmark && params.iteration[] > 5

    if params.discharge_voltage_IE > 0
        ie_index = argmin(abs.(grid.cell_centers .- params.IE_position))
        # compute sheath drops on each face (plasma minus metal potential)

        # metal potential of the IE (electrode)
        V_IE_wall = V_L - params.discharge_voltage_IE # e.g. discharge 500 V, discharge_voltage_IE 100 -> IE metal at 400 V

        # plasma-side potentials adjacent to the IE on left and right
        use_ie_sheath = false
        if use_ie_sheath
            Vs_L = ie_sheath_potential(cache, ie_index, :L)
            Vs_R = ie_sheath_potential(cache, ie_index, :R)
            V_IEp_L = V_IE_wall + Vs_L
            V_IEp_R = V_IE_wall - Vs_R
        else
            Vs_L = NaN
            Vs_R = NaN
            V_IEp_L = V_IE_wall
            V_IEp_R = V_IE_wall
        end
        @printf("  Sheath Potential: Vs: %.3f Vs_L: %.3f Vs_R: %.3f\n", Vs[], Vs_L, Vs_R)
    else
        ie_index = 0
        V_IE_wall = NaN
        V_IEp_L = NaN
        V_IEp_R = NaN
        Vs_L = NaN
        Vs_R = NaN
    end

    if ie_index != 0
        # left branch: integrate from anode (V_L) to left-side plasma potential at IE (V_IEp_L)
        cell_range_A = 1:ie_index
        Id_L_IE[] = integrate_discharge_current(grid, cache, V_L, V_IEp_L, "left", apply_drag, cell_range_A)
        
        # right branch: integrate from right-side plasma potential at IE (V_IEp_R) to cathode (V_R)
        cell_range_B = ie_index:length(grid.cell_centers)
        Id_IE_R[] = integrate_discharge_current(grid, cache, V_IEp_R, V_R, "right", apply_drag, cell_range_B)
        
        Id[] = NaN
    else
        Id_L_IE[] = NaN
        Id_IE_R[] = NaN
        cell_range = 1:length(grid.cell_centers) # use the full domain if no IE is present
        Id[] = integrate_discharge_current(grid, cache, V_L, V_R, "both", apply_drag, cell_range)
    end
    @printf("  Discharge current: Id: %.3f Id_L_IE: %.3f Id_IE_R: %.3f\n", Id[], Id_L_IE[], Id_IE_R[])
    @printf("  ie_index: %d\n", ie_index)

    # Compute electric field and potential
    update_electric_field!(∇ϕ, ie_index, cache, apply_drag) # ∇ϕ is updated in place from electric field solve
    
    # Integrate potential, pass the plasma-side potential on the right face of the IE
    integrate_potential!(ϕ, ∇ϕ, grid, V_L, V_IE_wall, ie_index) # integrates ϕ from ∇ϕ and applies boundary conditions, updates ϕ in place

    # Compute the electron velocity and electron kinetic energy
    @inbounds for i in eachindex(ue)
        # je + ji = Id / A
        if i > ie_index && params.discharge_voltage_IE > 0
            ue[i] = (ji[i] - Id_IE_R[] / channel_area[i]) / e / ne[i]
        elseif i <= ie_index && params.discharge_voltage_IE > 0
            ue[i] = (ji[i] - Id_L_IE[] / channel_area[i]) / e / ne[i]
        else
            ue[i] = (ji[i] - Id[] / channel_area[i]) / e / ne[i]
        end
    end

    # Kinetic energy in both axial and azimuthal directions is accounted for
    electron_kinetic_energy!(K, νe, B, ue)
    return
end

# Compute the axially-constant discharge current using Ohm's law
# Eq. 19 from "Numerical and Experimental Investigation of Longitudinal Oscillations in Hall Thrusters"
# https://www.mdpi.com/2226-4310/8/6/148
function integrate_discharge_current(grid, cache, V_L, V_R, edgeSide, apply_drag, cell_range)
    (; ∇pe, μ, ne, ji, channel_area, avg_neutral_vel, avg_ion_vel, νei, νen, νan) = cache

    # context
    # ΔV = ∫ [ji/(e·ne·μ) + ∇pe/(e·ne)] dz  -  Id · ∫ [1/(e·ne·μ·A)] dz
    # integrand_1 = ji/(e·μ·ne) + ∇pe/(e·ne)
    # integrant_2 = 1/(e·ne·μ·A)
    # ΔV = int1 - Id * int2,
    # Id = (ΔV + int1) / int2, solve for Id since it's the unknown

    # Compute integrands at all cell centers
    integrand_1 = cache.cell_cache_1
    integrand_2 = cache.cell_cache_2

    # assuming indices correspond to cell centers
    # only compute over specified cell range, to allow separate integration over IE and non-IE regions
    @inbounds for i in cell_range
        integrand_1[i] = (ji[i] / e / μ[i] + ∇pe[i]) / ne[i]
        integrand_2[i] = inv(e * ne[i] * μ[i] * channel_area[i])

        if (apply_drag)
            ion_drag_1 = avg_ion_vel[i] * (νei[i] + νan[i]) * me / e
            neutral_drag_1 = avg_neutral_vel[i] * νen[i] * me / e
            integrand_1[i] -= ion_drag_1 + neutral_drag_1
        end
    end

    # Replace values at edges of the cell range being integrated over
    if (edgeSide == "left" || edgeSide == "both")
        integrand_1[first(cell_range)] = 0.5 * (integrand_1[first(cell_range)] + integrand_1[first(cell_range)+1])
        integrand_2[first(cell_range)] = 0.5 * (integrand_2[first(cell_range)] + integrand_2[first(cell_range)+1])
    end
    if (edgeSide == "right" || edgeSide == "both")
        integrand_1[last(cell_range)]  = 0.5 * (integrand_1[last(cell_range)-1] + integrand_1[last(cell_range)])
        integrand_2[last(cell_range)]  = 0.5 * (integrand_2[last(cell_range)-1] + integrand_2[last(cell_range)])
    end

    # Compute integrals using trapezoidal rule around edges
    int1 = 0.0
    int2 = 0.0
    @inbounds for (idx, i) in enumerate(cell_range[1:end-1])

        # initially assume edge values are the same as the adjacent cell centers
        zL = grid.cell_centers[i]
        zR = grid.cell_centers[i + 1]

        f1_L = integrand_1[i]
        f1_R = integrand_1[i + 1]

        f2_L = integrand_2[i]
        f2_R = integrand_2[i + 1]

        # adjust z value at the edges
        if i == first(cell_range) && (edgeSide == "left" || edgeSide == "both")
            # if first cell in the range on the left edge, use the left z value
            zL = grid.edges[1]
        end
        if i == last(cell_range)-1 && (edgeSide == "right" || edgeSide == "both")
            # if last cell in the range on the right edge, use the right z value
            zR = grid.edges[end]
        end

        int1 += 0.5 * (zR - zL) * (f1_L + f1_R)
        int2 += 0.5 * (zR - zL) * (f2_L + f2_R)
    end

    ΔV = V_L - V_R

    I = (ΔV + int1) / int2

    return I
end

function electron_kinetic_energy!(K, νe, B, ue)
    # K = 1/2 me ue^2
    #   = 1/2 me (ue^2 + ue_θ^2)
    #   = 1/2 me (ue^2 + Ωe^2 ue^2)
    #   = 1/2 me (1 + Ωe^2) ue^2
    #   divide by e to get units of eV
    return @. K = 0.5 * me * (1 + (e * B / me / νe)^2) * ue^2 / e
end
