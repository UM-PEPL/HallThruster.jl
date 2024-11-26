# update useful quantities relevant for potential, electron energy and fluid solve
function update_electrons!(params, config, t = 0)
    (; grid, cache, simulation) = params
    (; nn, B, ue, Tev, ∇ϕ, ϕ, pe, ne, nϵ, μ, ∇pe, νan, νc, νen, νei, radial_loss_frequency,
    Z_eff, νiz, νex, νe, ji, Id, νew_momentum, κ, Vs, K,
    anom_multiplier, channel_area) = cache

    L_ch = config.thruster.geometry.channel_length

    # Update electron temperature and pressure given new density
    update_temperature!(Tev, nϵ, ne, params.min_Te)
    update_pressure!(pe, nϵ, config.LANDMARK)

    # Update electron-ion, electron-neutral, and total classical collisions
    if (config.electron_ion_collisions)
        freq_electron_ion!(νei, ne, Tev, Z_eff)
    end

    freq_electron_neutral!(νen, params.electron_neutral_collisions, nn, Tev)
    freq_electron_classical!(νc, νen, νei, νiz, νex, config.LANDMARK)

    # Update wall collisions
    @inbounds for i in eachindex(radial_loss_frequency)
        # Compute wall collision frequency, with transition function to force no momentum wall collisions in plume
        radial_loss_frequency[i] = freq_electron_wall(
            config.wall_loss_model, params, i,)
        νew_momentum[i] = radial_loss_frequency[i] * linear_transition(
            grid.cell_centers[i], L_ch, config.transition_length, 1.0, 0.0,)
    end

    # Update anomalous transport
    t > 0 && config.anom_model(νan, params, config)

    # Smooth anomalous transport model
    if config.anom_smoothing_iters > 0
        smooth!(νan, params.cache.cell_cache_1, iters = config.anom_smoothing_iters)
    end

    @inbounds for i in eachindex(νan)
        # Multiply by anom anom multiplier for PID control
        νan[i] *= anom_multiplier[]

        # Compute total collision frequency and electron mobility
        νe[i] = νc[i] + νan[i] + νew_momentum[i]
        μ[i] = electron_mobility(νe[i], B[i])
    end

    # Compute anode sheath potential
    Vs[] = anode_sheath_potential(params)

    # Compute the discharge current by integrating the momentum equation over the whole domain
    V_L = config.discharge_voltage + Vs[]
    V_R = config.cathode_potential
    un = config.neutral_velocity
    apply_drag = !config.LANDMARK && params.iteration[] > 5
    Id[] = integrate_discharge_current(grid, cache, V_L, V_R, un, apply_drag)

    # Compute the electron velocity and electron kinetic energy
    @inbounds for i in eachindex(ue)
        # je + ji = Id / A
        ue[i] = (ji[i] - Id[] / channel_area[i]) / e / ne[i]
    end

    # Kinetic energy in both axial and azimuthal directions is accounted for
    electron_kinetic_energy!(K, νe, B, ue)

    # Compute potential gradient and pressure gradient
    compute_pressure_gradient!(∇pe, pe, grid.cell_centers)

    # Compute electric field
    compute_electric_field!(∇ϕ, cache, un, apply_drag)

    # update electrostatic potential and potential gradient on edges
    cumtrapz!(ϕ, grid.cell_centers, ∇ϕ, config.discharge_voltage + Vs[])

    #update thermal conductivity
    config.conductivity_model(κ, params)

    # Update the electron temperature and pressure
    update_electron_energy!(params, config, params.dt[])

    # Update the anomalous collision frequency multiplier to match target current
    anom_multiplier[] = exp(apply_controller(
        simulation.current_control, Id[], log(anom_multiplier[]), params.dt[],
    ))
end

# Compute the axially-constant discharge current using Ohm's law
function integrate_discharge_current(grid, cache, V_L, V_R, un, apply_drag)
    (; ∇pe, μ, ne, ji, channel_area) = cache

    ncells = length(grid.cell_centers)

    int1 = 0.0
    int2 = 0.0

    if (apply_drag)
        (; νei, νen, νan, ui) = cache
    end

    @inbounds for i in 1:(ncells - 1)
        Δz = grid.dz_edge[i]

        int1_1 = (ji[i] / e / μ[i] + ∇pe[i]) / ne[i]
        int1_2 = (ji[i + 1] / e / μ[i + 1] + ∇pe[i + 1]) / ne[i + 1]

        if (apply_drag)
            ion_drag_1 = ui[1, i] * (νei[i] + νan[i]) * me / e
            ion_drag_2 = ui[1, i + 1] * (νei[i + 1] + νan[i + 1]) * me / e
            neutral_drag_1 = un * νen[i] * me / e
            neutral_drag_2 = un * νen[i + 1] * me / e
            int1_1 -= ion_drag_1 + neutral_drag_1
            int1_2 -= ion_drag_2 + neutral_drag_2
        end

        int1 += 0.5 * Δz * (int1_1 + int1_2)

        int2_1 = inv(e * ne[i] * μ[i] * channel_area[i])
        int2_2 = inv(e * ne[i + 1] * μ[i + 1] * channel_area[i + 1])

        int2 += 0.5 * Δz * (int2_1 + int2_2)
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
    @. K = 0.5 * me * (1 + (e * B / me / νe)^2) * ue^2 / e
end

function compute_pressure_gradient!(∇pe, pe, z_cell)
    ncells = length(z_cell)

    # Pressure gradient (forward)
    ∇pe[1] = forward_difference(pe[1], pe[2], pe[3], z_cell[1], z_cell[2], z_cell[3])

    # Centered difference in interior cells
    @inbounds for j in 2:(ncells - 1)
        # Compute pressure gradient
        ∇pe[j] = central_difference(
            pe[j - 1], pe[j], pe[j + 1], z_cell[j - 1], z_cell[j], z_cell[j + 1],)
    end

    # pressure gradient (backward)
    ∇pe[end] = backward_difference(
        pe[end - 2], pe[end - 1], pe[end], z_cell[end - 2], z_cell[end - 1], z_cell[end],)

    return nothing
end
