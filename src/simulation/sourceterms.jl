function apply_reactions!(fluid_arr, params)
    (;
        ionization_reactions,
        ionization_reactant_indices,
        ionization_product_indices,
        cache, landmark,
    ) = params

    rxns = zip(
        ionization_reactions, ionization_reactant_indices, ionization_product_indices,
    )

    return apply_reactions!(fluid_arr, rxns, cache, landmark)
end

function apply_reactions!(fluids, rxns, cache, landmark)
    (; inelastic_losses, νiz, ϵ, ne, K) = cache

    # Zero ionization frequency and inelastic losses and compute electron density
    @inbounds begin
        # Update electron density (TODO check if this is optimal)
        ne .= 0.0
        for fluid in fluids
            for i in eachindex(ne)
                ne[i] += fluid.species.Z * fluid.density[i] / fluid.species.element.m
            end
        end
        νiz .= 0.0
        inelastic_losses .= 0.0
        @. ϵ = cache.nϵ / cache.ne
        if !landmark
            @. ϵ += K
        end
    end

    dt_max = Inf
    for (rxn, reactant_index, product_index) in rxns
        reactant = fluids[reactant_index]
        product = fluids[product_index]
        _dt = apply_reaction!(reactant, product, ne, ϵ, rxn, νiz, inelastic_losses, landmark)
        dt_max = min(_dt, dt_max)
    end

    cache.dt_iz[] = dt_max
    return
end

function apply_reaction!(reactant, product, ne, ϵ, rxn, νiz, inelastic_losses, landmark)
    dt_max = Inf
    reactant_velocity = reactant.const_velocity
    inv_m = 1 / reactant.species.element.m

    @inbounds @simd for i in 2:(length(reactant.density) - 1)
        r = rate_coeff(rxn, ϵ[i])
        ρ_reactant = reactant.density[i]
        ρdot = reaction_rate(r, ne[i], ρ_reactant)
        dt_max = min(dt_max, ρ_reactant / ρdot)
        ndot = ρdot * inv_m
        νiz[i] += ndot / ne[i]

        inelastic_losses[i] += ndot * rxn.energy

        # Change in density due to ionization
        reactant.dens_ddt[i] -= ρdot
        product.dens_ddt[i] += ρdot

        if !landmark
            # Momentum transfer due to ionization
            if reactant.type != _ContinuityOnly
                reactant_velocity = reactant.momentum[i] / ρ_reactant
                reactant.mom_ddt[i] -= ρdot * reactant_velocity
            end

            product.mom_ddt[i] += ρdot * reactant_velocity
        end
    end

    return dt_max
end

@inline reaction_rate(rate_coeff, ne, n_reactant) = rate_coeff * ne * n_reactant

function apply_ion_acceleration!(fluids::Vector{FluidContainer}, grid, cache)
    dt_max = Inf

    @inbounds for fluid in fluids
        Z = fluid.species.Z
        m = fluid.species.element.m
        qe_m = Z * e / m

        @simd for i in 2:(length(fluid.dens_ddt) - 1)
            qE_m = -qe_m * cache.∇ϕ[i]
            dz = grid.dz_cell[i]

            Q_accel = qE_m * fluid.density[i]
            dt_max = min(dt_max, abs(dz / qE_m))
            fluid.mom_ddt[i] += Q_accel
        end
    end

    return cache.dt_E[] = sqrt(dt_max)
end

function apply_ion_wall_losses!(fluid_containers, params)
    (; thruster, cache, grid, transition_length, wall_loss_scale) = params
    (; continuity, isothermal) = fluid_containers

    geometry = thruster.geometry
    L_ch = geometry.channel_length
    inv_Δr = inv(geometry.outer_radius - geometry.inner_radius)
    h = wall_loss_scale * edge_to_center_density_ratio()

    neutral_fluid = continuity[1]
    return @inbounds for ion_fluid in isothermal
        m = ion_fluid.species.element.m
        Z = ion_fluid.species.Z
        qe_m = Z * e / m

        for i in 2:(length(ion_fluid.density) - 1)
            u_bohm = sqrt(qe_m * cache.Tev[i])
            in_channel = linear_transition(grid.cell_centers[i], L_ch, transition_length, 1.0, 0.0)
            νiw = in_channel * u_bohm * inv_Δr * h

            density_loss = ion_fluid.density[i] * νiw
            momentum_loss = ion_fluid.momentum[i] * νiw

            # Neutrals gain density due to ion recombination at the walls
            neutral_fluid.dens_ddt[i] += density_loss
            ion_fluid.dens_ddt[i] -= density_loss
            ion_fluid.mom_ddt[i] -= momentum_loss
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
