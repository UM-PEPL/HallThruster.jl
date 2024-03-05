function apply_reactions!(dU::AbstractArray{T}, U::AbstractArray{T}, params, i::Int64) where T
    (;index, ionization_reactions, index, ionization_reactant_indices, ionization_product_indices, cache) = params
    (;inelastic_losses, νiz) = cache

    ne = electron_density(U, params, i)
    K = if params.config.LANDMARK
        0.0
    else
        params.cache.K[i]
    end

    ϵ  = cache.nϵ[i] / ne + K

    dt_max = Inf
    νiz[i] = 0.0
    inelastic_losses[i] = 0.0

    @inbounds for (rxn, reactant_inds, product_inds) in zip(ionization_reactions, ionization_reactant_indices, ionization_product_indices)
        product_index = product_inds[]
        rate_coeff = rxn.rate_coeff(ϵ)

        for reactant_index in reactant_inds
            ρ_reactant = U[reactant_index, i]
            ρdot = reaction_rate(rate_coeff, ne, ρ_reactant)
            ndot = ρdot / params.config.propellant.m
            νiz[i] += ndot / ne
            inelastic_losses[i] += ndot * rxn.energy
            dt_max = min(dt_max, ρ_reactant / ρdot)

            # Change in density due to ionization
            dU[reactant_index, i] -= ρdot
            dU[product_index, i]  += ρdot

            if !params.config.LANDMARK
                # Momentum transfer due to ionization
                if reactant_index == index.ρn[1]
                    reactant_velocity = params.config.neutral_velocity
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
    E = -cache.∇ϕ[i]
    mi = params.config.propellant.m
    dt_max = Inf

    Δx = z_edge[right_edge(i)] - z_edge[left_edge(i)]

    @inbounds for Z in 1:config.ncharge
        Q_accel = Z * e * U[index.ρi[Z], i] / mi * E
        dt_max = min(dt_max, sqrt(mi * Δx / Z / e / abs(E)))
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

function excitation_losses!(U, params, i)
    (;excitation_reactions, cache, excitation_reactant_indices) = params
    (;νex, Tev, inelastic_losses) = cache
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

    νex[i] = 0.0
    @inbounds for (reactant_inds, rxn) in zip(excitation_reactant_indices, excitation_reactions)
        rate_coeff = rxn.rate_coeff(ϵ)
        for reactant_index in reactant_inds
            n_reactant = U[reactant_index, i] / mi
            ndot = reaction_rate(rate_coeff, ne, n_reactant)
            inelastic_losses[i] += ndot * (rxn.energy - K)
            νex[i] += ndot / ne
        end
    end

    return inelastic_losses[i]
end

function electron_kinetic_energy(U, params, i)
    (;ue, Ωe²) = params.cache
    return 0.5 *  me * (1 + Ωe²[i]) * ue[i]^2 / e
end

function source_electron_energy(U, params, i)
    (;ne, ue, ∇ϕ) = params.cache

    # Compute ohmic heating term, which is the rate at which energy is transferred out of the electron
    # drift (kinetic energy) into thermal energy
    ohmic_heating = ne[i] * ue[i] * ∇ϕ[i]

    # add excitation losses to total inelastic losses
    excitation_losses!(U, params, i)

    # compute wall losses
    wall_loss = ne[i] * wall_power_loss(params.config.wall_loss_model, U, params, i)

    params.cache.wall_losses[i] = wall_loss
    params.cache.ohmic_heating[i] = ohmic_heating

    return ohmic_heating - wall_loss - params.cache.inelastic_losses[i]
end
