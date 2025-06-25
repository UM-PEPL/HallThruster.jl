function apply_reactions!(dU, U, params)
    (;
        index, ionization_reactions, index,
        ionization_reactant_indices, ionization_product_indices,
        cache, mi, landmark, ncharge, neutral_velocity,
    ) = params

    rxns = zip(
        ionization_reactions, ionization_reactant_indices, ionization_product_indices,
    )

    return apply_reactions!(dU, U, cache, index, ncharge, mi, landmark, neutral_velocity, rxns)
end

function apply_reactions!(dU, U, cache, index, ncharge, mi, landmark, un, rxns)
    (; inelastic_losses, νiz, ϵ, ne, K) = cache

    inv_m = inv(mi)
    ncells = length(ne)

    νiz .= 0.0
    inelastic_losses .= 0.0
    @inbounds for i in eachindex(ne)
        ne[i] = 0.0
        for Z in 1:ncharge
            ne[i] += Z * U[index.ρi[Z], i]
        end
        ne[i] *= inv_m
    end
    inv_ne = cache.cell_cache_1
    @. inv_ne = inv(cache.ne)
    @. ϵ = cache.nϵ * inv_ne
    if !landmark
        @. ϵ += K
    end

    dt_max = Inf

    @inbounds for (rxn, reactant_index, product_index) in rxns
        for i in 2:(ncells - 1)
            r = rate_coeff(rxn, ϵ[i])
            ρ_reactant = U[reactant_index, i]
            ρdot = reaction_rate(r, ne[i], ρ_reactant)
            dt_max = min(dt_max, ρ_reactant / ρdot)
            ndot = ρdot * inv_m
            νiz[i] += ndot * inv_ne[i]

            inelastic_losses[i] += ndot * rxn.energy

            # Change in density due to ionization
            dU[reactant_index, i] -= ρdot
            dU[product_index, i] += ρdot

            if !landmark
                # Momentum transfer due to ionization
                if reactant_index == index.ρn
                    reactant_velocity = un
                else
                    reactant_velocity = U[reactant_index + 1, i] / ρ_reactant
                    dU[reactant_index + 1, i] -= ρdot * reactant_velocity
                end

                dU[product_index + 1, i] += ρdot * reactant_velocity
            end
        end
    end

    return cache.dt_iz[] = dt_max
end

@inline reaction_rate(rate_coeff, ne, n_reactant) = rate_coeff * ne * n_reactant

function apply_ion_acceleration!(dU, U, params)
    (; index, grid, cache, mi, ncharge) = params
    return apply_ion_acceleration!(dU, U, grid, cache, index, mi, ncharge)
end
function apply_ion_acceleration!(dU, U, grid, cache, index, mi, ncharge)
    inv_m = inv(mi)
    inv_e = inv(e)
    dt_max = Inf

    @inbounds for i in 2:(length(grid.cell_centers) - 1)
        E = -cache.∇ϕ[i]
        Δz = grid.dz_cell[i]
        inv_E = inv(abs(E))

        @inbounds for Z in 1:(ncharge)
            Q_accel = Z * e * U[index.ρi[Z], i] * inv_m * E
            dt_max = min(dt_max, sqrt(mi * Δz * inv_e * inv_E / Z))
            dU[index.ρiui[Z], i] += Q_accel
        end
    end

    return cache.dt_E[] = dt_max
end

function apply_ion_wall_losses!(dU, U, params)
    (; index, ncharge, mi, thruster, cache, grid, transition_length, wall_loss_scale) = params
    geometry = thruster.geometry
    L_ch = geometry.channel_length
    inv_Δr = inv(geometry.outer_radius - geometry.inner_radius)
    e_inv_m = e / mi

    h = wall_loss_scale * edge_to_center_density_ratio()

    return @inbounds for i in 2:(length(grid.cell_centers) - 1)
        u_bohm = sqrt(e_inv_m * cache.Tev[i])
        in_channel = linear_transition(
            grid.cell_centers[i], L_ch, transition_length, 1.0, 0.0,
        )
        νiw_base = in_channel * u_bohm * inv_Δr * h

        for Z in 1:ncharge
            νiw = sqrt(Z) * νiw_base
            density_loss = U[index.ρi[Z], i] * νiw
            momentum_loss = U[index.ρiui[Z], i] * νiw

            # Neutrals gain density due to ion recombination at the walls
            dU[index.ρi[Z], i] -= density_loss
            dU[index.ρn, i] += density_loss
            dU[index.ρiui[Z], i] -= momentum_loss
        end
    end
end

function apply_user_ion_source_terms!(
        dU, U, params, source_neutrals, source_ion_continuity, source_ion_momentum,
    )
    (; index, grid, ncharge) = params
    # Apply user-provided source terms
    ncells = length(grid.cell_centers)
    return @inbounds for i in 2:(ncells - 1)
        # User-provided neutral source term
        dU[index.ρn, i] += source_neutrals(U, params, i)
        for Z in 1:ncharge
            # User-provided ion source terms
            dU[index.ρi[Z], i] += source_ion_continuity[Z](U, params, i)
            dU[index.ρiui[Z], i] += source_ion_momentum[Z](U, params, i)
        end
    end
end

function excitation_losses!(Q, cache, landmark, grid, excitation_reactions)
    (; νex, ϵ, nn, ne, K) = cache
    ncells = length(grid.cell_centers)

    @. νex = 0.0
    @inbounds for rxn in excitation_reactions
        for i in 2:(ncells - 1)
            r = rate_coeff(rxn, ϵ[i])
            ndot = reaction_rate(r, ne[i], nn[i])
            νex[i] += ndot / ne[i]
            Q[i] += ndot * (rxn.energy - !landmark * K[i])
        end
    end

    return nothing
end

function ohmic_heating!(Q, cache, landmark)
    (; ne, ue, ∇ϕ, K, νe, ue, ∇pe) = cache
    # Compute ohmic heating term, which is the rate at which energy is transferred out of the electron
    # drift (kinetic energy) into thermal energy
    if (landmark)
        # Neglect kinetic energy, so the rate of energy transfer into thermal energy is equivalent to
        # the total input power into the electrons (j⃗ ⋅ E⃗ = -mₑnₑ|uₑ|²νₑ)
        @. Q = ne * ue * ∇ϕ
    else
        # Do not neglect kinetic energy, so ohmic heating term is mₑnₑ|uₑ|²νₑ + ue ∇pe = 2nₑKνₑ + ue ∇pe
        # where K is the electron bulk kinetic energy, 1/2 * mₑ|uₑ|²
        @. Q = 2 * ne * K * νe + ue * ∇pe
    end
    return nothing
end

function source_electron_energy!(Q, params, wall_loss_model)
    (; cache, landmark, grid, excitation_reactions) = params
    (; ne, ohmic_heating, wall_losses, inelastic_losses) = cache

    # compute ohmic heating
    ohmic_heating!(ohmic_heating, cache, landmark)

    # add excitation losses to total inelastic losses
    excitation_losses!(
        inelastic_losses, cache, landmark, grid, excitation_reactions,
    )

    # compute wall losses
    wall_power_loss!(wall_losses, wall_loss_model, params)

    # Compute net energy source, i.e heating minus losses
    @. Q = ohmic_heating - ne * wall_losses - inelastic_losses

    return Q
end
