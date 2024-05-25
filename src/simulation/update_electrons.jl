# update useful quantities relevant for potential, electron energy and fluid solve
function update_electrons!(U, params, t = 0)
    (;index, control_current, target_current, Kp, Ti, mi, ncells) = params
    (;
        B, ue, Tev, ∇ϕ, ϕ, pe, ne, nϵ, μ, ∇pe, νan, νc, νen, νei, radial_loss_frequency,
        Z_eff, νiz, νex, νe, ji, Id, νew_momentum, κ, ni, ui, Vs, nn, niui,
        Id_smoothed, smoothing_time_constant, anom_multiplier,
        errors, channel_area
    ) = params.cache

    # Allow for system interrupts
    yield()

    # Update the current iteration
    params.iteration[1] += 1

    # Apply fluid boundary conditions
    @views left_boundary_state!(U[:, 1], U, params)
    @views right_boundary_state!(U[:, end], U, params)

    # Update plasma quantities
    @inbounds for i in 1:ncells
        # Compute number density for each neutral fluid
        nn[i] = U[index.ρn, i] / params.config.propellant.m

        # Compute ion derived quantities
        ne[i] = 0.0
        Z_eff[i] = 0.0
        ji[i] = 0.0
        @inbounds for Z in 1:params.config.ncharge
            _ni = U[index.ρi[Z], i] / mi
            _niui = U[index.ρiui[Z], i] / mi
            ni[Z, i] = _ni
            niui[Z, i] = _niui
            ui[Z, i] = _niui / _ni
            ne[i] += Z * _ni
            Z_eff[i] += _ni
            ji[i] += Z * e * _niui
        end

        # Compute electron number density, making sure it is above floor
        ne[i] = max(params.config.min_number_density, ne[i])

        # Effective ion charge state (density-weighted average charge state)
        Z_eff[i] = max(1.0, ne[i] / Z_eff[i])

        # Compute new electron temperature
        Tev[i] = 2/3 * max(params.config.min_electron_temperature, nϵ[i]/ne[i])

        # Compute electron pressure
        pe[i] = if params.config.LANDMARK
            # The LANDMARK benchmark uses nϵ instead of pe in the potential solver, but we use pe, so
            # we need to define pe = 3/2 ne Tev
            3/2 * ne[i] * Tev[i]
        else
            # Otherwise, just use typical ideal gas law.
            ne[i] * Tev[i]
        end
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
        νen[i] = freq_electron_neutral(params.electron_neutral_collisions, nn[i], Tev[i])

        # Compute total classical collision frequency
        # If we're not running the LANDMARK benchmark, include momentum transfer due to inelastic collisions
        νc[i] = νen[i] + νei[i] + !params.config.LANDMARK * (νiz[i] + νex[i])

        # Compute wall collision frequency, with transition function to force no momentum wall collisions in plume
        radial_loss_frequency[i] = freq_electron_wall(params.config.wall_loss_model, params, i)
        νew_momentum[i] =  radial_loss_frequency[i]* params.config.transition_function(params.z_cell[i], params.L_ch, 1.0, 0.0)

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
    Id[] = _discharge_current(params)

    # Compute the electron velocity and electron kinetic energy
    @inbounds for i in 1:ncells
        # je + ji = Id / A
        ue[i] = (ji[i] - Id[] / channel_area[i]) / e / ne[i]

        # Kinetic energy in both axial and azimuthal directions is accounted for
        params.cache.K[i] = electron_kinetic_energy(params, i)
    end

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

    # Update the anomalous collision frequency multiplier to match target
    # discharge current
    if control_current && t > 0
        Ki = Kp / Ti

        A1 = Kp + Ki*dt
        A2 = -Kp

        α = 1 - exp(-params.dt[]/smoothing_time_constant[])
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


# TODO: differentiate this from the postprocessing one
function _discharge_current(params)
    (;cache, Δz_edge, ϕ_L, ϕ_R, ncells) = params
    (;∇pe, μ, ne, ji, Vs, channel_area) = cache

    int1 = 0.0
    int2 = 0.0

    @inbounds for i in 1:ncells-1
        Δz = Δz_edge[i]

        int1_1 = (ji[i] / e / μ[i] + ∇pe[i]) / ne[i]
        int1_2 = (ji[i+1] / e / μ[i+1] + ∇pe[i+1]) / ne[i+1]

        int1 += 0.5 * Δz * (int1_1 + int1_2)

        int2_1 = inv(e * ne[i] * μ[i] * channel_area[i])
        int2_2 = inv(e * ne[i+1] * μ[i+1] * channel_area[i+1])

        int2 += 0.5 * Δz * (int2_1 + int2_2)
    end

    ΔV = ϕ_L - ϕ_R + Vs[]

    I = (ΔV + int1) / int2

    return I
end

function compute_electric_field!(∇ϕ, params)
    (;cache) = params
    (;ji, Id, ne, μ, ∇pe, channel_area) = cache

    for i in eachindex(∇ϕ)
        ∇ϕ[i] = -((Id[] / channel_area[i] - ji[i]) / e / μ[i] - ∇pe[i]) / ne[i]
    end

    return ∇ϕ
end

function compute_pressure_gradient!(∇pe, params)
    (; pe) = params.cache
    (;z_cell, ncells) = params

    # Pressure gradient (forward)
    ∇pe[1] = forward_difference(pe[1], pe[2], pe[3], z_cell[1], z_cell[2], z_cell[3])

    # Centered difference in interior cells
    @inbounds for j in 2:ncells-1
        # Compute pressure gradient
        ∇pe[j] = central_difference(pe[j-1], pe[j], pe[j+1], z_cell[j-1], z_cell[j], z_cell[j+1])
    end

    # pressure gradient (backward)
    ∇pe[end] = backward_difference(pe[end-2], pe[end-1], pe[end], z_cell[end-2], z_cell[end-1], z_cell[end])

    return nothing
end

function smooth!(x, x_cache; iters=1)
    if iters > 0
        x_cache .= x
        x_cache[1] = x[2]
        x_cache[end-1] = x[end]
        for i in 2:length(x)-1
            if i == 2 || i == length(x)-1
                x_cache[i] = 0.5 * x[i] + 0.25 * (x[i-1] + x[i+1])
            else
                x_cache[i] = 0.4 * x[i] + 0.24 * (x[i-1] + x[i+1]) + 0.06 * (x[i-2] + x[i+2])
            end
        end
        x .= x_cache
        smooth!(x, x_cache; iters = iters-1)
    else
        return x
    end
end
