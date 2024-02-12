
abstract type InitialCondition end

struct DefaultInitialization <: InitialCondition end

initialize!(U, params) = initialize!(U, params, params.config.initial_condition)

function initialize!(U, params, model::InitialCondition)
    throw(ArgumentError("Function HallThruster.initialize!(U, params, model::$(typeof(model)) not yet implemented. For InitialCondition types other than DefaultInitialization(), this must be defined by the user!"))
end

function initialize!(U, params, ::DefaultInitialization)
    (;z_cell, config, index) = params
    (;ncharge, anode_Te, cathode_Te, domain, thruster, propellant, discharge_voltage, anode_mass_flow_rate) = config
    mi = propellant.m
    L_ch = thruster.geometry.channel_length
    z0 = domain[1]

    ni_center = L_ch / 2
    ni_width = L_ch / 3
    ni_min = 2e17
    ni_max = 1e18
    scaling_factor = sqrt(discharge_voltage / 300) * (anode_mass_flow_rate / 5e-6)
    ion_density_function = (z, Z) -> mi * scaling_factor * (ni_min + (ni_max - ni_min) * exp(-(((z-z0) - ni_center) / ni_width)^2)) / Z^2

    bohm_velocity = Z -> -sqrt(Z * e * anode_Te / mi)

    final_velocity = Z -> sqrt(2 * Z * e * discharge_voltage / mi)
    scale(Z) = 2/3 * (final_velocity(Z) - bohm_velocity(Z))
    ion_velocity_f1(z, Z) = bohm_velocity(Z) + scale(Z) * ((z-z0) / L_ch)^2
    ion_velocity_f2(z, Z) = lerp(z, z0 + L_ch, domain[2], ion_velocity_f1(L_ch, Z), final_velocity(Z))

    ion_velocity_function = (z, Z) -> if z - z0 < L_ch
        ion_velocity_f1(z, Z)
    else
        ion_velocity_f2(z, Z)
    end

    ρn_0 = inlet_neutral_density(config)
    # add recombined neutrals
    for Z in 1:config.ncharge
        ρn_0 -= ion_velocity_function(0.0, Z) * ion_density_function(0.0, Z) / config.neutral_velocity
    end

    # Beam neutral density at outlet
    ρn_1 = 0.01 * ρn_0
    neutral_function = z -> SmoothIf(transition_length = L_ch / 6)(z-z0, L_ch / 2, ρn_0, ρn_1)

    number_density_function = z -> sum(Z * ion_density_function(z, Z) / mi for Z in 1:ncharge)

    Te_baseline = z -> lerp(z, domain[1], domain[2], anode_Te, cathode_Te)
    Te_min = min(anode_Te, cathode_Te)
    Te_max = (config.discharge_voltage / 10)
    Te_width = L_ch/3

    energy_function = z -> 3/2 * (Te_baseline(z) + (Te_max - Te_min) * exp(-(((z-z0) - L_ch) / Te_width)^2))

    # Fill the state vector
    for (i, z) in enumerate(z_cell)

        U[index.ρn[1], i] = neutral_function(z)

        for Z in 1:params.config.ncharge
            U[index.ρi[Z], i] = ion_density_function(z, Z)
            U[index.ρiui[Z], i] = ion_density_function(z, Z) * ion_velocity_function(z, Z)
        end

        U[index.nϵ, i] = number_density_function(z) * energy_function(z)
    end

    return U
end

function initialize_dt!(U, params)
    
    
    """
    Flux dt 
    """
    (;index, config, cache, fluids, num_neutral_fluids, z_edge, Δz_edge, 
        ionization_reactions, ionization_reactant_indices, CFL
    ) = params
    (;scheme, LANDMARK, propellant, electron_pressure_coupled, ncharge, 
        min_electron_temperature, min_number_density
    ) = config
    (;
        B, ue, Tev, ∇ϕ, ϕ, pe, ne, μ, ∇pe, νan, νc, νen, νei, νew,
        Z_eff, νiz, νex, νe, ji, Id, νew, ni, ui, Vs, nn, nn_tot, niui, K,
        anom_multiplier, channel_area, dt_u, dt_E, dt_iz, UL, UR
    ) = cache

    #constants
    mi = propellant.m
    if LANDMARK
        Te_fac = 1.0
    else
        Te_fac = 2/3
    end

    #mesh properties
    coupled = electron_pressure_coupled
    (nvars,  ncells) = size(U)
    nedges = ncells-1


    #first need to compute left and right edge states
    @inbounds for i in 2:ncells-1
        for j in 1:nvars
            if scheme.reconstruct

                is_velocity_index = false
                for Z in 1:params.config.ncharge
                    if j == index.ρiui[Z]
                        is_velocity_index = true
                        break
                    end
                end

                if is_velocity_index # reconstruct velocity as primitive variable instead of momentum density
                    u₋ = U[j, i-1]/U[j-1, i-1]
                    uᵢ = U[j, i]/U[j-1, i]
                    u₊ = U[j, i+1]/U[j-1, i+1]
                    uR, uL = reconstruct(u₋, uᵢ, u₊, scheme.limiter)

                    ρL = UL[j-1, right_edge(i)] #use previously-reconstructed edge density to compute momentum
                    ρR = UR[j-1, left_edge(i)]
                    UL[j, right_edge(i)] = uL*ρL
                    UR[j, left_edge(i)] = uR*ρR
                else
                    u₋ = U[j, i-1]
                    uᵢ = U[j, i]
                    u₊ = U[j, i+1]

                    UR[j, left_edge(i)], UL[j, right_edge(i)] = reconstruct(u₋, uᵢ, u₊, scheme.limiter)
                end
            else
                UL[j, right_edge(i)] = U[j, i]
                UR[j, left_edge(i)]  = U[j, i]
            end
        end
    end

    #apply boundary conditions is default true, so just use the true version here
    @views left_boundary_state!(UL[:, 1], U, params)
    @views right_boundary_state!(UR[:, end], U, params)


    @inbounds for i in 1:nedges

        # Compute number density
        neL = 0.0
        neR = 0.0
        for Z in 1:ncharge
            niL = UL[index.ρi[Z], i] / mi
            niR = UR[index.ρi[Z], i] / mi
            neL = neL + Z * niL
            neR = neR + Z * niR
        end

        Δz = Δz_edge[i]

        # Compute electron temperature
        TeL = max(min_electron_temperature, Te_fac * UL[index.nϵ, i] / neL)
        TeR = max(min_electron_temperature, Te_fac * UR[index.nϵ, i] / neR)

        for Z in 1:ncharge
            fluid_ind = Z + num_neutral_fluids
            fluid = fluids[fluid_ind]
            γ = fluid.species.element.γ
            UL_ions = (UL[index.ρi[Z], i], UL[index.ρiui[Z], i],)
            UR_ions = (UR[index.ρi[Z], i], UR[index.ρiui[Z], i],)

            uL = velocity(UL_ions, fluid)
            TL = temperature(UL_ions, fluid)
            uR = velocity(UR_ions, fluid)
            TR = temperature(UR_ions, fluid)

            # If we're using the electron-pressure-coupled method, then the sound speed includes a contribution from
            # the electron temperature and approaches the ion acoustic speed
            charge_factor = Z * e * coupled

            # Sound speeds
            aL = sqrt((charge_factor * TeL + γ * kB * TL) / mi)
            aR = sqrt((charge_factor * TeR + γ * kB * TR) / mi)

            # There are several possible waves, with speeds, u-a, u+a, and u, both left- and right-running
            λ₁₂⁺ = 0.5 * (uR + abs(uR))
            λ₁₂⁻ = 0.5 * (uL - abs(uL))
            λ₃⁺  = 0.5 * (uR + aR  + abs(uR + aR))
            λ₃⁻  = 0.5 * (uL + aL - abs(uL + aL))
            λ₄⁺  = 0.5 * (uR - aR  + abs(uR - aR))
            λ₄⁻  = 0.5 * (uL - aL  - abs(uL - aL))

            # Maximum of all of these wave speeds
            s_max = max(abs(λ₁₂⁺), abs(λ₁₂⁻), abs(λ₃⁺), abs(λ₃⁻), abs(λ₄⁺), abs(λ₄⁻))

            # a Δt / Δx = 1 for CFL condition, user-supplied CFL number restriction applied later, in update_values
            dt_max = Δz / s_max
            dt_u[i] = dt_max
        end
    end



    """
    Acceleration dt 
    """

    # loop over cells 
    @inbounds for i in 1:ncells

        # Compute neutral number densities for each neutral fluid
        nn_tot[i] = 0.0
        for j in 1:num_neutral_fluids
            nn[j, i] = U[index.ρn[j], i] / mi
            nn_tot[i] += nn[j, i]
        end

        # Compute ion densities and velocities
        for Z in 1:ncharge
            ni[Z, i] = U[index.ρi[Z], i] / mi
            ui[Z, i] = U[index.ρiui[Z], i] / U[index.ρi[Z], i]
            niui[Z, i] = U[index.ρiui[Z], i] / mi
        end

        # Compute electron number density, making sure it is above floor
        ne[i] = max(min_number_density, electron_density(U, params, i))

        # Same with electron temperature
        Tev[i] = 2/3 * max(min_electron_temperature, U[index.nϵ, i]/ne[i])

        pe[i] = if LANDMARK
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
        if !LANDMARK
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
        K[i] = electron_kinetic_energy(U, params, i)
    end

    # Compute potential gradient and pressure gradient
    compute_pressure_gradient!(∇pe, params)

    # Compute electric field
    compute_electric_field!(∇ϕ, params)

    #perform the computation of the limit based on the potential gradient 
    @inbounds for i in 2:ncells-1
        dt_max = Inf
        Δx = z_edge[right_edge(i)] - z_edge[left_edge(i)]

        @inbounds for Z in 1:ncharge

            dt_max = min(dt_max, sqrt(mi * Δx / Z / e / abs(∇ϕ[i])))

        end

        dt_E[i] = dt_max
    end
   
    """
    Ionization dt
    """

    @inbounds for i in 2:ncells-1
        K_cell = if LANDMARK
            0.0
        else
            K[i]
        end
    
        ϵ  = U[index.nϵ, i] / ne[i] + K_cell
    
        dt_max = Inf
    
        @inbounds for (rxn, reactant_inds) in zip(ionization_reactions, ionization_reactant_indices)
            
            rate_coeff = rxn.rate_coeff(ϵ)
    
            for reactant_index in reactant_inds
                ρ_reactant = U[reactant_index, i]
    
                ρdot = reaction_rate(rate_coeff, ne[i], ρ_reactant)
    
                dt_max = min(dt_max, ρ_reactant / ρdot)
    
            end
        end
    
        dt_iz[i] = dt_max
    end

    """
    Total dt
    """
    @inbounds for i in 2:ncells-1
        left = left_edge(i)
        right = right_edge(i)
        
        params.cache.dt_cell[i] = min(
                sqrt(CFL) * params.cache.dt_E[i],
                CFL * params.cache.dt_iz[i],
                CFL * params.cache.dt_u[left],
                CFL * params.cache.dt_u[right]
            )

        params.cache.dt[] = min(params.cache.dt_cell[i], params.cache.dt[])
    end
end