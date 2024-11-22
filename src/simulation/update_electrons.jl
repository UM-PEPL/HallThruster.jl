
# update useful quantities relevant for potential, electron energy and fluid solve
function update_electrons!(params, t = 0)
    (; config, grid, cache, simulation) = params
    (; B, ue, Tev, ∇ϕ, ϕ, pe, ne, nϵ, μ, ∇pe, νan, νc, νen, νei, radial_loss_frequency,
    Z_eff, νiz, νex, νe, ji, Id, νew_momentum, κ, Vs, K,
    anom_multiplier, channel_area) = cache

    ncells = length(grid.cell_centers)
    L_ch = config.thruster.geometry.channel_length

    # Update plasma quantities based on new density
    pe_factor = config.LANDMARK ? 1.5 : 1.0
    @inbounds for i in 1:ncells
        # Compute new electron temperature
        Tev[i] = 2 / 3 * max(params.min_Te, nϵ[i] / ne[i])
        # Compute electron pressure
        pe[i] = pe_factor * ne[i] * Tev[i]
    end

    # Update electron-ion collisions
    if (params.config.electron_ion_collisions)
        @inbounds for i in 1:ncells
            νei[i] = freq_electron_ion(ne[i], Tev[i], Z_eff[i])
        end
    end

    # Update other collisions
    @inbounds for i in 1:ncells
        # Compute electron-neutral and electron-ion collision frequencies
        νen[i] = freq_electron_neutral(params, i)

        # Compute total classical collision frequency
        # If we're not running the LANDMARK benchmark, include momentum transfer due to inelastic collisions
        νc[i] = νen[i] + νei[i] + !config.LANDMARK * (νiz[i] + νex[i])

        # Compute wall collision frequency, with transition function to force no momentum wall collisions in plume
        radial_loss_frequency[i] = freq_electron_wall(
            params.config.wall_loss_model, params, i,)
        νew_momentum[i] = radial_loss_frequency[i] * linear_transition(
            grid.cell_centers[i], L_ch, params.config.transition_length, 1.0, 0.0,)
    end

    # Update anomalous transport
    t > 0 && params.config.anom_model(νan, params)

    # Smooth anomalous transport model
    if params.config.anom_smoothing_iters > 0
        smooth!(νan, params.cache.cell_cache_1, iters = params.config.anom_smoothing_iters)
    end

    @inbounds for i in 1:ncells
        # Multiply by anom anom multiplier for PID control
        νan[i] *= anom_multiplier[]

        # Compute total collision frequency and electron mobility
        νe[i] = νc[i] + νan[i] + νew_momentum[i]
        μ[i] = electron_mobility(νe[i], B[i])
    end

    # Compute anode sheath potential
    Vs[] = anode_sheath_potential(params)

    # Compute the discharge current by integrating the momentum equation over the whole domain
    Id[] = integrate_discharge_current(params)

    # Compute the electron velocity and electron kinetic energy
    @inbounds for i in 1:ncells
        # je + ji = Id / A
        ue[i] = (ji[i] - Id[] / channel_area[i]) / e / ne[i]
    end

    # Kinetic energy in both axial and azimuthal directions is accounted for
    electron_kinetic_energy!(K, params)

    # Compute potential gradient and pressure gradient
    compute_pressure_gradient!(∇pe, params)

    # Compute electric field
    compute_electric_field!(∇ϕ, params)

    # update electrostatic potential and potential gradient on edges
    solve_potential!(ϕ, params)

    #update thermal conductivity
    params.config.conductivity_model(κ, params)

    # Update the electron temperature and pressure
    update_electron_energy!(params, params.dt[])

    # Update the anomalous collision frequency multiplier to match target current
    anom_multiplier[] = exp(apply_controller(
        simulation.current_controller, Id[], log(anom_multiplier[]), params.dt[],
    ))
end

# Compute the axially-constant discharge current using Ohm's law
function integrate_discharge_current(params)
    (; grid, config, cache, iteration) = params
    (; ∇pe, μ, ne, ji, Vs, channel_area) = cache

    V_L, V_R = config.discharge_voltage, config.cathode_potential
    ncells = length(grid.cell_centers)

    int1 = 0.0
    int2 = 0.0

    apply_drag = false & !config.LANDMARK & (iteration[] > 5)

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
            neutral_drag_1 = params.config.neutral_velocity * νen[i] * me / e
            neutral_drag_2 = params.config.neutral_velocity * νen[i + 1] * me / e
            int1_1 -= ion_drag_1 + neutral_drag_1
            int1_2 -= ion_drag_2 + neutral_drag_2
        end

        int1 += 0.5 * Δz * (int1_1 + int1_2)

        int2_1 = inv(e * ne[i] * μ[i] * channel_area[i])
        int2_2 = inv(e * ne[i + 1] * μ[i + 1] * channel_area[i + 1])

        int2 += 0.5 * Δz * (int2_1 + int2_2)
    end

    ΔV = V_L - V_R + Vs[]

    I = (ΔV + int1) / int2

    return I
end

function compute_electric_field!(∇ϕ, params)
    (; config, cache, iteration) = params
    (; ji, Id, ne, μ, ∇pe, channel_area, ui, νei, νen, νan) = cache

    apply_drag = false & !config.LANDMARK & (iteration[] > 5)

    if (apply_drag)
        (; νei, νen, νan, ui) = cache
    end

    for i in eachindex(∇ϕ)
        E = ((Id[] / channel_area[i] - ji[i]) / e / μ[i] - ∇pe[i]) / ne[i]

        if (apply_drag)
            ion_drag = ui[1, i] * (νei[i] + νan[i]) * me / e
            neutral_drag = config.neutral_velocity * (νen[i]) * me / e
            E += ion_drag + neutral_drag
        end

        ∇ϕ[i] = -E
    end

    return ∇ϕ
end

function electron_kinetic_energy!(K, params)
    (; νe, B, ue) = params.cache
    # K = 1/2 me ue^2
    #   = 1/2 me (ue^2 + ue_θ^2)
    #   = 1/2 me (ue^2 + Ωe^2 ue^2)
    #   = 1/2 me (1 + Ωe^2) ue^2
    #   divide by e to get units of eV
    @. K = 0.5 * me * (1 + (e * B / me / νe)^2) * ue^2 / e
end

function compute_pressure_gradient!(∇pe, params)
    (; pe) = params.cache
    (; grid) = params
    z_cell = grid.cell_centers
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

function smooth!(x, x_cache; iters = 1)
    if iters > 0
        x_cache .= x
        x_cache[1] = x[2]
        x_cache[end - 1] = x[end]
        for i in 2:(length(x) - 1)
            if i == 2 || i == length(x) - 1
                x_cache[i] = 0.5 * x[i] + 0.25 * (x[i - 1] + x[i + 1])
            else
                x_cache[i] = 0.4 * x[i] + 0.24 * (x[i - 1] + x[i + 1]) +
                             0.06 * (x[i - 2] + x[i + 2])
            end
        end
        x .= x_cache
        smooth!(x, x_cache; iters = iters - 1)
    else
        return x
    end
end
