# update useful quantities relevant for potential, electron energy and fluid solve
function update_electrons!(params, config, t = 0)
    (; Tev, ne, nϵ, νan, νc, νen, νei, radial_loss_frequency, Z_eff, νiz, νex, νew_momentum, κ) = params.cache
    (; source_energy, wall_loss_model, conductivity_model, anom_model) = config

    # Update electron temperature given new density
    update_temperature!(Tev, nϵ, ne, params.min_Te)

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
    update_walls!(radial_loss_frequency, νew_momentum, wall_loss_model, params)

    # Update anomalous transport
    t > 0 && anom_model(νan, params)

    # Update mobility, discharge current, potential, and electron velocity
    update_electrical_vars!(params)

    # update the thermal conductivity
    conductivity_model(κ, params)

    # Update the electron energy density,  temperature and pressure
    update_electron_energy!(params, wall_loss_model, source_energy, params.dt[])

    return
end

function update_walls!(radial_loss_frequency, νew_momentum, wall_loss_model, params)
    (; thruster, grid, transition_length) = params
    L_ch = thruster.geometry.channel_length

    freq_electron_wall!(radial_loss_frequency, wall_loss_model, params)

    # Update wall collisions
    @inbounds for i in eachindex(radial_loss_frequency)
        # Compute wall collision frequency, with transition function to force no momentum wall collisions in plume
        νew_momentum[i] = radial_loss_frequency[i] * linear_transition(
            grid.cell_centers[i], L_ch, transition_length, 1.0, 0.0,
        )
    end
    return nothing
end

function update_electrical_vars!(params)
    (; cache, anom_smoothing_iters, landmark, grid) = params
    (;
        cell_cache_1, νan, νe, νc, μ, B, νew_momentum, anom_multiplier,
        Vs, ue, ji, channel_area, ne, Id, Vd, K,  ϕ, ∇ϕ,
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
    V_L = Vd[] + Vs[]
    V_R = params.cathode_coupling_voltage

    apply_drag = !landmark && params.iteration[] > 5

    Id[] = integrate_discharge_current(grid, cache, V_L, V_R, apply_drag)
    Vd[] = update_circuit(params.filter_circuit, params.discharge_voltage, Id[], params.dt[])

    # Compute electric field and potential
    update_electric_field!(∇ϕ, cache, apply_drag)
    integrate_potential!(ϕ, ∇ϕ, grid, V_L)

    # Compute the electron velocity and electron kinetic energy
    @inbounds for i in eachindex(ue)
        # je + ji = Id / A
        ue[i] = (ji[i] - Id[] / channel_area[i]) / e / ne[i]
    end

    # Kinetic energy in both axial and azimuthal directions is accounted for
    electron_kinetic_energy!(K, νe, B, ue)
    return
end

# Compute the axially-constant discharge current using Ohm's law
function integrate_discharge_current(grid, cache, V_L, V_R, apply_drag)
    (; ∇pe, μ, ne, ji, channel_area, avg_neutral_vel, avg_ion_vel, νei, νen, νan) = cache

    # Compute integrands at all cell centers
    integrand_1 = cache.cell_cache_1
    integrand_2 = cache.cell_cache_2

    @inbounds for i in eachindex(grid.cell_centers)
        integrand_1[i] = (ji[i] / e / μ[i] + ∇pe[i]) / ne[i]
        integrand_2[i] = inv(e * ne[i] * μ[i] * channel_area[i])

        if (apply_drag)
            ion_drag_1 = avg_ion_vel[i] * (νei[i] + νan[i]) * me / e
            neutral_drag_1 = avg_neutral_vel[i] * νen[i] * me / e
            integrand_1[i] -= ion_drag_1 + neutral_drag_1
        end
    end

    # Replace left and right values with edge values
    integrand_1[1] = 0.5 * (integrand_1[1] + integrand_1[2])
    integrand_1[end] = 0.5 * (integrand_1[end - 1] + integrand_1[end])
    integrand_2[1] = 0.5 * (integrand_2[1] + integrand_2[2])
    integrand_2[end] = 0.5 * (integrand_2[end - 1] + integrand_2[end])

    # Compute integrals using trapezoidal rule around edges
    int1 = 0.0
    int2 = 0.0
    @inbounds for (i, z_edge) in enumerate(grid.edges)
        zL = grid.cell_centers[i]
        zR = grid.cell_centers[i + 1]

        f1_L = integrand_1[i]
        f1_R = integrand_1[i + 1]

        f2_L = integrand_2[i]
        f2_R = integrand_2[i + 1]

        # account for boundary cells
        if i == 1 || i == length(grid.edges)
            zL = z_edge
        elseif i == length(grid.edges)
            zR = z_edge
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
