function apply_reactions!(dU::AbstractArray{T}, U::AbstractArray{T}, params, i::Int64) where T
    (;index, ionization_reactions, index, ionization_reactant_indices, ionization_product_indices) = params

    ne = electron_density(U, params, i)
    K = if params.config.LANDMARK
        0.0
    else
        params.cache.K[i]
    end

    ϵ  = U[index.nϵ, i] / ne + K

    dt_max = Inf

    @inbounds for (rxn, reactant_inds, product_inds) in zip(ionization_reactions, ionization_reactant_indices, ionization_product_indices)
        product_index = product_inds[]

        rate_coeff = rxn.rate_coeff(ϵ)

        for reactant_index in reactant_inds
            ρ_reactant = U[reactant_index, i]

            ρdot = reaction_rate(rate_coeff, ne, ρ_reactant)

            dt_max = min(dt_max, ρ_reactant / ρdot)

            # Change in density due to ionization
            dU[reactant_index, i] -= ρdot
            dU[product_index, i]  += ρdot

            if !params.config.LANDMARK
                # Momentum transfer due to ionization
                if reactant_index == index.ρn[1]
                    reactant_velocity = params.config.neutral_velocity
                elseif reactant_index == index.ρn[2]
                    reactant_velocity = params.background_neutral_velocity
                else
                    reactant_velocity = U[reactant_index + 1, i] / U[reactant_index, i]
                    dU[reactant_index + 1, i] -= ρdot * reactant_velocity
                end

                dU[product_index + 1, i] += ρdot * reactant_velocity
            end
        end
    end

    params.cache.dt_iz[i] = dt_max
end

@inline reaction_rate(rate_coeff, ne, n_reactant) = rate_coeff * ne * n_reactant
@inline reaction_rate(rxn, ne, n_reactant, ϵ) = rxn.rate_coeff(ϵ) * n_reactant * ne

function apply_ion_acceleration!(dU, U, params, i)
    (;cache, config, index, z_edge) = params
    (;∇ϕ, ue, μ) = cache
    coupled = config.electron_pressure_coupled
    mi = params.config.propellant.m

    dt_max = Inf

    Δx = z_edge[right_edge(i)] - z_edge[left_edge(i)]

    @inbounds for Z in 1:config.ncharge

        Q_coupled = ue[i] / μ[i]
        Q_uncoupled = ∇ϕ[i]

        Q_accel = -Z * e * U[index.ρi[Z], i] / mi * (coupled * Q_coupled + (1 - coupled) * Q_uncoupled)

        dt_max = min(dt_max, sqrt(mi * Δx / Z / e / abs(∇ϕ[i])))

        dU[index.ρiui[Z], i] += Q_accel
    end

    params.cache.dt_E[i] = dt_max
end

function apply_ion_wall_losses!(dU, U, params, i)
    (;index, config, z_cell, L_ch) = params
    (;ncharge, propellant, wall_loss_model, thruster) = config

    geometry = thruster.geometry
    Δr = geometry.outer_radius - geometry.inner_radius

    mi = propellant.m

    if wall_loss_model isa WallSheath
        α = wall_loss_model.α
    elseif wall_loss_model isa NoWallLosses
        return
    else
        α = 1.0
    end

    @inbounds for Z in 1:ncharge
        
        in_channel = params.config.transition_function(z_cell[i], L_ch, 1.0, 0.0)
        u_bohm = sqrt(Z * e * params.cache.Tev[i] / mi)
        νiw = α * in_channel * u_bohm / Δr

        density_loss  = U[index.ρi[Z], i] * νiw
        momentum_loss = U[index.ρiui[Z], i] * νiw
        
        dU[index.ρi[Z],   i] -= density_loss
        dU[index.ρiui[Z], i] -= momentum_loss
        
        # Neutrals gain density due to ion recombination at the walls
        dU[index.ρn[1], i] += density_loss
    end
end

function source_electron_energy!(Q, U, params, i)
    Q[params.index.nϵ] = source_electron_energy(U, params, i)
end

function inelastic_losses!(U, params, i)
    (;ionization_reactions, excitation_reactions, cache, ionization_reactant_indices, excitation_reactant_indices) = params
    (;νex, νiz, Tev) = cache
    inelastic_loss = 0.0
    mi = params.config.propellant.m
    ne = cache.ne[i]

    K = if params.config.LANDMARK
        # Neglect electron kinetic energy if LANDMARK
        0.0
    else
        # Include kinetic energy contribution
        cache.K[i]
    end

    # Total electron energy
    ϵ = 3/2 * Tev[i] + K

    νiz[i] = 0.0
    νex[i] = 0.0

    @inbounds for (reactant_inds, rxn) in zip(ionization_reactant_indices, ionization_reactions)
        rate_coeff = rxn.rate_coeff(ϵ)
        for reactant_index in reactant_inds
            n_reactant = U[reactant_index, i] / mi
            ndot = reaction_rate(rate_coeff, ne, n_reactant)
            inelastic_loss += ndot * rxn.energy
            νiz[i] += ndot / ne
        end
    end

    @inbounds for (reactant_inds, rxn) in zip(excitation_reactant_indices, excitation_reactions)
        rate_coeff = rxn.rate_coeff(ϵ)
        for reactant_index in reactant_inds
            n_reactant = U[reactant_index, i] / mi
            ndot = reaction_rate(rate_coeff, ne, n_reactant)
            inelastic_loss += ndot * (rxn.energy - K)
            νex[i] += ndot / ne
        end
    end

    return inelastic_loss
end

function electron_kinetic_energy(U, params, i)
    νe = params.cache.νe[i]
    B  = params.cache.B[i]
    ue = params.cache.ue[i]
    Ωe = e * B / me / νe
    return 0.5 *  me * (1 + Ωe^2) * ue^2 / e
end

function source_electron_energy(U, params, i)
    ne = params.cache.ne[i]
    ue = params.cache.ue[i]
    ∇ϕ = params.cache.∇ϕ[i]

    # Compute ohmic heating term, which is the rate at which energy is transferred out of the electron
    # drift (kinetic energy) into thermal inergy
    if params.config.LANDMARK
        # Neglect kinetic energy, so the rate of energy transfer into thermal energy is equivalent to
        # the total input power into the electrons (j⃗ ⋅ E⃗ = -mₑnₑ|uₑ|²νₑ)
        ohmic_heating  = ne * ue * ∇ϕ
    else
        # Do not neglect kinetic energy, so ohmic heating term is mₑnₑ|uₑ|²νₑ + ue ∇pe = 2nₑKνₑ + ue ∇pe
        # where K is the electron bulk kinetic energy, 1/2 * mₑ|uₑ|²
        K = params.cache.K[i]
        νe = params.cache.νe[i]

        z = params.z_cell
        pe = params.cache.pe

        # Upwind difference for pressure gradient
        if ue > 0
            ∇pe = (pe[i] - pe[i-1]) / (z[i] - z[i-1])
        else
            ∇pe = (pe[i+1] - pe[i]) / (z[i+1] - z[i])
        end

        ohmic_heating = 2 * ne * K * νe + ue * ∇pe
    end

    wall_loss      = ne * wall_power_loss(params.config.wall_loss_model, U, params, i)
    inelastic_loss = inelastic_losses!(U, params, i)

    params.cache.wall_losses[i] = wall_loss
    params.cache.inelastic_losses[i] = inelastic_loss
    params.cache.ohmic_heating[i] = ohmic_heating

    return ohmic_heating - wall_loss - inelastic_loss
end
