# update useful quantities relevant for potential, electron energy and fluid solve
function update_electrons!(U, params, t = 0)
    (;index, control_current, target_current, Kp, Ti) = params
    (;
        B, ue, Tev, ∇ϕ, ϕ, pe, ne, μ, ∇pe, νan, νc, νen, νei, νew,
        Z_eff, νiz, νex, νe, ji, Id, νew, ni, ui, Vs, nn, nn_tot, niui,
        Id_smoothed, smoothing_time_constant, anom_multiplier,
        errors, channel_area
    ) = params.cache

    # Update the current iteration
    params.iteration[1] += 1

    # Apply fluid boundary conditions
    @views left_boundary_state!(U[:, 1], U, params)
    @views right_boundary_state!(U[:, end], U, params)

    ncells = size(U, 2)

    t_sub = 0.0
    num_subcycles = 0

    while t_sub < params.dt[]
        num_subcycles += 1

        # Update electron quantities
        @inbounds for i in 1:ncells

            # Compute neutral number densities for each neutral fluid
            nn_tot[i] = 0.0
            for j in 1:params.num_neutral_fluids
                nn[j, i] = U[index.ρn[j], i] / params.config.propellant.m
                nn_tot[i] += nn[j, i]
            end

            # Compute ion densities and velocities
            for Z in 1:params.config.ncharge
                ni[Z, i] = U[index.ρi[Z], i] / params.config.propellant.m
                ui[Z, i] = U[index.ρiui[Z], i] / U[index.ρi[Z], i]
                niui[Z, i] = U[index.ρiui[Z], i] / params.config.propellant.m
            end

            # Compute electron number density, making sure it is above floor
            ne[i] = max(params.config.min_number_density, electron_density(U, params, i))

            # Same with electron temperature
            Tev[i] = 2/3 * max(params.config.min_electron_temperature, U[index.nϵ, i]/ne[i])

            pe[i] = if params.config.LANDMARK
                # The LANDMARK benchmark uses nϵ instead of pe in the potential solver, but we use pe, so
                # we need to define pe = 3/2 ne Tev
                3/2 * ne[i] * Tev[i]
            else
                # Otherwise, just use typical ideal gas law.
                ne[i] * Tev[i]
            end
            # Compute electron-neutral and electron-ion collision frequencies
            νen[i] = freq_electron_neutral(U, params, i)
            νei[i] = freq_electron_ion(U, params, i)

            # Compute total classical collision frequency
            νc[i] = νen[i] + νei[i]
            if !params.config.LANDMARK
                # Add momentum transfer due to ionization and excitation
                νc[i] += νiz[i] + νex[i]
            end

            # Compute anomalous collision frequency and wall collision frequencies
            νew[i] = freq_electron_wall(params.config.wall_loss_model, U, params, i)
        end

        # Update anomalous transport
        params.config.anom_model(νan, params)

        # Smooth anomalous transport model
        if params.config.anom_smoothing_iters > 0
            smooth!(νan, params.cache.cell_cache_1, iters = params.config.anom_smoothing_iters)
        end

        @inbounds for i in 1:ncells
            # Multiply by anom anom multiplier for PID control
            νan[i] *= anom_multiplier[]

            # Compute total collision frequency and electron mobility
            νe[i] = νc[i] + νan[i] + νew[i]
            μ[i] = electron_mobility(νe[i], B[i])

            # Effective ion charge state (density-weighted average charge state)
            Z_eff[i] = compute_Z_eff(U, params, i)

            # Ion current
            ji[i] = ion_current_density(U, params, i)
        end

        # Compute anode sheath potential
        Vs[] = anode_sheath_potential(U, params)

        # Compute the discharge current by integrating the momentum equation over the whole domain
        Id[] = discharge_current(U, params)

        # Compute the electron velocity and electron kinetic energy
        @inbounds for i in 1:ncells
            # je + ji = Id / A
            ue[i] = (ji[i] - Id[] / channel_area[i]) / e / ne[i]

            # Kinetic energy in both axial and azimuthal directions is accounted for
            params.cache.K[i] = electron_kinetic_energy(U, params, i)
        end

        # Compute potential gradient and pressure gradient
        compute_pressure_gradient!(∇pe, params)

        # Compute electric field
        compute_electric_field!(∇ϕ, params)

        # update electrostatic potential and potential gradient on edges
        solve_potential_cell!(ϕ, params)

        dt_min = Inf

        @inbounds for i in 2:ncells-1
            Q = source_electron_energy(U, params, i)
            dt_min = min(dt_min, abs(params.CFL * 3 * ne[i] * Tev[i] / Q))
        end

        dt_sub = min(dt_min, params.dt[] - t_sub)

        # Update the electron temperature and pressure
        update_electron_energy!(U, params, dt_sub)

        t_sub += dt_sub
    end

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


function discharge_current(U::Array, params)
    (;cache, Δz_edge, ϕ_L, ϕ_R) = params
    (;∇pe, μ, ne, ji, Vs, channel_area) = cache

    ncells = size(U, 2)

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
    (;z_cell) = params

    ncells = length(z_cell)

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
