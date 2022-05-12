function apply_reactions!(dU::AbstractArray{T}, U::AbstractArray{T}, params, i::Int64) where T
    (;index, ionization_reactions, index, ionization_reactant_indices, ionization_product_indices) = params

    ne = electron_density(U, params, i)
    K = if params.config.LANDMARK
        0.0
    else
        params.cache.K[i]
    end

    ϵ  = U[index.nϵ, i] / ne + K

    @inbounds for (rxn, reactant_index, product_index) in zip(ionization_reactions, ionization_reactant_indices, ionization_product_indices)
        ρ_reactant = U[reactant_index, i]

        ρdot = reaction_rate(rxn, ne, ρ_reactant, ϵ)

        # Change in density due to ionization
        dU[reactant_index, i] -= ρdot
        dU[product_index, i]  += ρdot

        if !params.config.LANDMARK
            # Momentum transfer due to ionization
            if reactant_index != index.ρn
                reactant_velocity = U[reactant_index + 1, i] / U[reactant_index, i]
                dU[reactant_index + 1, i] -= ρdot * reactant_velocity
            else
                reactant_velocity = params.config.neutral_velocity
            end
            dU[product_index + 1, i] += ρdot * reactant_velocity
        end

    end
end

@inline reaction_rate(rxn, ne, n_reactant, ϵ) = rxn.rate_coeff(ϵ) * n_reactant * ne

function apply_ion_acceleration!(dU, U, params, i)
    (;cache, config, index) = params
    (;∇ϕ, ue, μ) = cache
    coupled = config.electron_pressure_coupled
    mi = params.config.propellant.m

    @inbounds for Z in 1:config.ncharge

        Q_coupled = ue[i] / μ[i]
        Q_uncoupled = ∇ϕ[i]

        Q_accel = -Z * e * U[index.ρi[Z], i] / mi * (coupled * Q_coupled + (1 - coupled) * Q_uncoupled)
        dU[index.ρiui[Z], i] += Q_accel
    end

end

function apply_plume_losses!(dU, U, params, i)
    (;index, config, z_cell, L_ch, fluids) = params
    (;ncharge, propellant, thruster) = config

    Δr = thruster.geometry.outer_radius - thruster.geometry.inner_radius
    mi = propellant.m

    if z_cell[i] > L_ch
        # Radial neutral losses
        Tn = temperature(SA[U[index.ρn, i]], fluids[1])
        dU[index.ρn, i] -= U[index.ρn, i] * sqrt(kB * Tn / 2 / π / mi) / Δr

        # Radial ion losses
        u_bohm = sqrt(e * params.cache.Tev[i] / mi)

        @inbounds for Z in 1:ncharge
            dU[index.ρi[Z],   i] -= U[index.ρi[Z],   i] * sqrt(Z) * u_bohm / Δr
            dU[index.ρiui[Z], i] -= U[index.ρiui[Z], i] * sqrt(Z) * u_bohm / Δr
        end
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

    @inbounds for (reactant_index, rxn) in zip(ionization_reactant_indices, ionization_reactions)
        n_reactant = U[reactant_index, i] / mi
        ndot = reaction_rate(rxn, ne, n_reactant, ϵ)
        inelastic_loss += ndot * rxn.energy
        νiz[i] += ndot / ne
    end

    @inbounds for (reactant_index, rxn) in zip(excitation_reactant_indices, excitation_reactions)
        n_reactant = U[reactant_index, i] / mi
        ndot = reaction_rate(rxn, ne, n_reactant, ϵ)
        inelastic_loss += ndot * (rxn.energy - K)
        νex[i] += ndot / ne
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

    wall_loss      = ne * params.config.wall_loss_model(U, params, i)
    inelastic_loss = inelastic_losses!(U, params, i)

    return ohmic_heating - wall_loss - inelastic_loss
end