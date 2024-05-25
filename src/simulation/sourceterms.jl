function apply_neutral_reaction!(dU, U, rxn, reactant_ind, product_ind, params)
    (;cache, config) = params
    (;ne, ϵ_tot, inelastic_losses, νiz) = params.cache

    # load some scratch space
    rate_coeff = cache.cell_cache_1
    n_dot = cache.cell_cache_2
    ρ_dot = cache.cell_cache_3
    ρu_dot = cache.cell_cache_4

    ρ_reactant =  @views U[reactant_ind, :]
    u_reactant = params.config.neutral_velocity
    dρ_reactant = @views dU[reactant_ind, :]
    dρ_product =  @views dU[product_ind, :]
    dρu_product = @views dU[product_ind + 1, :]

    # rate coefficient
    @. rate_coeff = rxn.rate_coeff(ϵ_tot)

    # change in reactant density
    @. ρ_dot = reaction_rate(rate_coeff, ne, ρ_reactant)
    @. n_dot = ρ_dot / config.propellant.m

    # change in reactant momentum
    @. ρu_dot = ρ_dot * u_reactant

    # apply changes
    @. dρ_reactant -= ρ_dot
    @. dρ_product  += ρ_dot
    @. dρu_product += ρu_dot

    # compute loss rates
    @. νiz += n_dot / ne
    @. inelastic_losses += n_dot * rxn.energy

    # update allowable timestep
    @. cache.dt_iz = min(cache.dt_iz, ρ_reactant / ρ_dot)
end

function apply_reaction!(dU, U, rxn, reactant_ind, product_ind, params)
    (;cache, config) = params
    (;ne, ϵ_tot, inelastic_losses, νiz) = params.cache # scratch space

    # load some scratch space
    rate_coeff = cache.cell_cache_1
    n_dot = cache.cell_cache_2
    ρ_dot = cache.cell_cache_3
    ρu_dot = cache.cell_cache_4

    ρ_reactant   = @views U[reactant_ind, :]
    ρu_reactant  = @views U[reactant_ind+1, :]
    dρ_reactant  = @views dU[reactant_ind, :]
    dρu_reactant = @views dU[reactant_ind + 1, :]
    dρ_product   = @views dU[product_ind, :]
    dρu_product  = @views dU[product_ind + 1, :]

    # rate coefficient
    @. rate_coeff = rxn.rate_coeff(ϵ_tot)

    # change in reactant density
    @. ρ_dot = reaction_rate(rate_coeff, ne, ρ_reactant)
    @. n_dot = ρ_dot / config.propellant.m

    # change in reactant momentum
    @. ρu_dot = ρ_dot * ρu_reactant / ρ_reactant

    # apply changes
    @. dρ_reactant -= ρ_dot
    @. dρ_product  += ρ_dot
    @. dρu_reactant -= ρu_dot
    @. dρu_product += ρu_dot

    # compute loss rates
    @. νiz += n_dot / ne
    @. inelastic_losses += n_dot * rxn.energy

    # update allowable timestep
    @. cache.dt_iz = min(cache.dt_iz, ρ_reactant / ρ_dot)
end

function apply_reactions!(dU::AbstractArray{T}, U::AbstractArray{T}, params) where T
    (;index, ionization_reactions, index, ionization_reactant_indices, ionization_product_indices, cache) = params

    # initialize some cache variables
    cache.νiz .= 0
    cache.inelastic_losses .= 0
    cache.dt_iz .= Inf

    @inbounds for (rxn, reactant_inds, product_inds) in zip(ionization_reactions, ionization_reactant_indices, ionization_product_indices)
        reactant_ind = reactant_inds[]
        product_ind  = product_inds[]
        if reactant_ind == params.index.ρn
            apply_neutral_reaction!(dU, U, rxn, reactant_ind, product_ind, params)
        else
            apply_reaction!(dU, U, rxn, reactant_ind, product_ind, params)
        end
    end
end

@inline reaction_rate(rate_coeff, ne, n_reactant) = rate_coeff * ne * n_reactant
@inline reaction_rate(rxn, ne, n_reactant, ϵ) = rxn.rate_coeff(ϵ) * n_reactant * ne

function apply_ion_acceleration!(dU, params)
    (;cache, config, index, ncells, Δz_cell) = params
    (;∇ϕ) = cache
    mi = params.config.propellant.m
    dt_max = Inf

    for i in 2:ncells-1
        E = ∇ϕ[i]
        @inbounds for Z in 1:config.ncharge
            Q_accel = -Z * e * cache.ni[Z, i] * E
            dt_max = min(dt_max, sqrt(mi * Δz_cell[i] / (Z * e * abs(E))))
            dU[index.ρiui[Z], i] += Q_accel
        end
        params.cache.dt_E[i] = dt_max
    end

    return nothing
end

function apply_ion_wall_losses!(dU, U, params)
    (;index, config, z_cell, L_ch, ncells) = params
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

    @inbounds for i in 2:ncells-1
        for Z in 1:ncharge
            in_channel = config.transition_function(z_cell[i], L_ch, 1.0, 0.0)
            u_bohm = sqrt(Z * e * params.cache.Tev[i] / mi)
            νiw = α * in_channel * u_bohm / Δr

            density_loss  = U[index.ρi[Z], i] * νiw
            momentum_loss = U[index.ρiui[Z], i] * νiw

            dU[index.ρi[Z],   i] -= density_loss
            dU[index.ρiui[Z], i] -= momentum_loss

            # Neutrals gain density due to ion recombination at the walls
            dU[index.ρn, i] += density_loss
        end
    end

    return nothing
end

function excitation_losses!(params, i)
    (;excitation_reactions, cache, excitation_reactant_indices) = params
    (;νex, Tev, inelastic_losses, nn) = cache
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
        for _ in reactant_inds
            n_reactant = nn[i]
            ndot = reaction_rate(rate_coeff, ne, n_reactant)
            inelastic_losses[i] += ndot * (rxn.energy - K)
            νex[i] += ndot / ne
        end
    end

    return inelastic_losses[i]
end

function electron_kinetic_energy(params, i)
    νe = params.cache.νe[i]
    B  = params.cache.B[i]
    ue = params.cache.ue[i]
    Ωe = e * B / me / νe
    return 0.5 *  me * (1 + Ωe^2) * ue^2 / e
end

function source_electron_energy(params, i)
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
        ∇pe = params.cache.∇pe[i]
        ohmic_heating = 2 * ne * K * νe + ue * ∇pe
    end

    # add excitation losses to total inelastic losses
    excitation_losses!(params, i)

    # compute wall losses
    wall_loss = ne * wall_power_loss(params.config.wall_loss_model, params, i)

    params.cache.wall_losses[i] = wall_loss
    params.cache.ohmic_heating[i] = ohmic_heating

    return ohmic_heating - wall_loss - params.cache.inelastic_losses[i]
end
