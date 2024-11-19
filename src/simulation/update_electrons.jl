# update useful quantities relevant for potential, electron energy and fluid solve
function update_electrons!(params, t = 0)
    (; control_current, target_current, Kp, Ti, pe_factor, ncells, cache, config) = params
    (; B, ue, Tev, electric_field, potential, pe, ne, nϵ, grad_pe,
    mobility, nu_e, nu_anom, nu_en, nu_ei, nu_class, nu_wall, nu_iz, nu_ex,
    radial_loss_frequency, Z_eff, ji, Id, κ, Vs, K,
    Id_smoothed, smoothing_time_constant, anom_multiplier,
    errors, channel_area) = cache

    # Update plasma quantities based on new density
    @inbounds for i in 1:ncells
        # Compute new electron temperature
        Tev[i] = 2 / 3 * max(params.Te_min, nϵ[i] / ne[i])
        # Compute electron pressure
        pe[i] = pe_factor * ne[i] * Tev[i]
    end

    # Update electron-ion collisions
    if (config.electron_ion_collisions)
        @inbounds for i in 1:ncells
            nu_ei[i] = freq_electron_ion(ne[i], Tev[i], Z_eff[i])
        end
    end

    # Update other collisions
    @inbounds for i in 1:ncells
        # Compute electron-neutral and electron-ion collision frequencies
        nu_en[i] = freq_electron_neutral(params, i)

        # Compute total classical collision frequency
        # If we're not running the LANDMARK benchmark, include momentum transfer due to inelastic collisions
        nu_class[i] = nu_en[i] + nu_ei[i] + !config.LANDMARK * (nu_iz[i] + nu_ex[i])

        # Compute wall collision frequency, with transition function to force no momentum wall collisions in plume
        radial_loss_frequency[i] = freq_electron_wall(
            config.wall_loss_model, params, i,)
        nu_wall[i] = radial_loss_frequency[i] * linear_transition(
            params.z_cell[i], params.L_ch, config.transition_length, 1.0, 0.0,)
    end

    # Update anomalous transport
    t > 0 && config.anom_model(nu_anom, params)

    # Smooth anomalous transport model
    if config.anom_smoothing_iters > 0
        smooth!(nu_anom, params.cache.cell_cache_1, iters = config.anom_smoothing_iters)
    end

    @inbounds for i in 1:ncells
        # Multiply by anom anom multiplier for PID control
        nu_anom[i] *= anom_multiplier[]

        # Compute total collision frequency and electron mobility
        nu_e[i] = nu_class[i] + nu_anom[i] + nu_wall[i]
        mobility[i] = electron_mobility(nu_e[i], B[i])
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
    compute_pressure_gradient!(grad_pe, params)

    # Compute electric field
    compute_electric_field!(electric_field, params)

    # update electrostatic potential and potential gradient on edges
    solve_potential!(potential, params)

    #update thermal conductivity
    config.conductivity_model(κ, params)

    # Update the electron temperature and pressure
    update_electron_energy!(params, params.dt[])

    # Update the anomalous collision frequency multiplier to match target
    # discharge current
    if control_current && t > 0
        Ki = Kp / Ti

        A1 = Kp + Ki * params.dt[]
        A2 = -Kp

        α = 1 - exp(-params.dt[] / smoothing_time_constant[])
        Id_smoothed[] = α * Id[] + (1 - α) * Id_smoothed[]

        errors[3] = errors[2]
        errors[2] = errors[1]
        errors[1] = target_current - Id_smoothed[]

        # PID controller
        if t > Ti
            new_anom_mult = log(anom_multiplier[]) + A1 * errors[1] + A2 * errors[2]
            anom_multiplier[] = exp(new_anom_mult)
        end
    elseif t == 0
        Id_smoothed[] = Id[]
        anom_multiplier[] = 1.0
    end
end

# Compute the axially-constant discharge current using Ohm's law
function integrate_discharge_current(params)
    (; cache, Δz_edge, V_L, V_R, ncells, iteration, config) = params
    (; grad_pe, mobility, ne, ji, Vs, channel_area) = cache

    int1 = 0.0
    int2 = 0.0

    apply_drag = false & !config.LANDMARK & (iteration[] > 5)

    if (apply_drag)
        (; νei, νen, νan, ui) = cache
    end
    un = config.neutral_velocity

    @inbounds for i in 1:(ncells - 1)
        Δz = Δz_edge[i]

        int1_1 = (ji[i] / e / mobility[i] + grad_pe[i]) / ne[i]
        int1_2 = (ji[i + 1] / e / mobility[i + 1] + grad_pe[i + 1]) / ne[i + 1]

        if (apply_drag)
            ion_drag_1 = ui[1, i] * (νei[i] + νan[i]) * me / e
            ion_drag_2 = ui[1, i + 1] * (νei[i + 1] + νan[i + 1]) * me / e
            neutral_drag_1 = un * νen[i] * me / e
            neutral_drag_2 = un * νen[i + 1] * me / e
            int1_1 -= ion_drag_1 + neutral_drag_1
            int1_2 -= ion_drag_2 + neutral_drag_2
        end

        int1 += 0.5 * Δz * (int1_1 + int1_2)

        int2_1 = inv(e * ne[i] * mobility[i] * channel_area[i])
        int2_2 = inv(e * ne[i + 1] * mobility[i + 1] * channel_area[i + 1])

        int2 += 0.5 * Δz * (int2_1 + int2_2)
    end

    ΔV = V_L - V_R + Vs[]

    I = (ΔV + int1) / int2

    return I
end

function compute_electric_field!(electric_field, params)
    (; cache, iteration, config) = params
    (; ji, Id, ne, mobility, grad_pe, channel_area, ui, nu_ei, nu_en, nu_anom) = cache

    apply_drag = false & !config.LANDMARK & (iteration[] > 5)

    un = config.neutral_velocity

    for i in eachindex(electric_field)
        E = ((Id[] / channel_area[i] - ji[i]) / e / mobility[i] - grad_pe[i]) / ne[i]

        if (apply_drag)
            ion_drag = ui[1, i] * (nu_ei[i] + nu_anom[i]) * me / e
            neutral_drag = un * (nu_en[i]) * me / e
            E += ion_drag + neutral_drag
        end

        electric_field[i] = E
    end

    return electric_field
end

function electron_kinetic_energy!(K, params)
    (; nu_e, B, ue) = params.cache
    # K = 1/2 me ue^2
    #   = 1/2 me (ue^2 + ue_θ^2)
    #   = 1/2 me (ue^2 + Ωe^2 ue^2)
    #   = 1/2 me (1 + Ωe^2) ue^2
    #   divide by e to get units of eV
    @. K = 0.5 * me * (1 + (e * B / me / nu_e)^2) * ue^2 / e
end

function compute_pressure_gradient!(∇pe, params)
    (; pe) = params.cache
    (; z_cell, ncells) = params

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
