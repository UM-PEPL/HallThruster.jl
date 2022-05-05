function apply_reactions!(dU::AbstractArray{T}, U::AbstractArray{T}, params, i::Int64) where T
    (;index, ionization_reactions, index, ionization_reactant_indices, ionization_product_indices) = params

    ne = electron_density(U, params, i)
    ϵ  = U[index.nϵ, i] / ne

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
    index = params.index
    (;∇ϕ, ue, μ) = params.cache
    coupled = params.config.electron_pressure_coupled
    mi = params.config.propellant.m

    @inbounds for Z in 1:params.config.ncharge
        # Ion acceleration
        Q_accel = coupled * ue[i] / μ[i] + (1 - coupled) * ∇ϕ[i]
        Q_accel = -Z * e * U[index.ρi[Z], i] / mi * Q_accel
        dU[index.ρiui[Z], i] += Q_accel

        # Radial ion loss
        #=
        geometry = params.config.thruster.geometry
        Δr = geometry.outer_radius - geometry.inner_radius
        u_bohm = sqrt(e * params.cache.Tev[i] / mi)
        dU[index.ρi[Z], i]   -= u_bohm / Δr * U[index.ρi[Z], i]
        dU[index.ρiui[Z], i] -= u_bohm / Δr * U[index.ρiui[Z], i]
        =#
    end
end

function source_electron_energy!(Q, U, params, i)
    Q[params.index.nϵ] = source_electron_energy(U, params, i)
end

function inelastic_losses!(U, params, i)
    (;ionization_reactions, excitation_reactions, cache, ionization_reactant_indices, excitation_reactant_indices) = params
    (;νex, νiz) = cache
    inelastic_loss = 0.0
    mi = params.config.propellant.m
    ne = cache.ne[i]
    ϵ  = 3/2 * cache.Tev[i]

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
        inelastic_loss += ndot * rxn.energy
        νex[i] += ndot / ne
    end

    return inelastic_loss
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
        # Do not neglect kinetic energy, so ohmic heating term is -mₑnₑ|uₑ|²νₑ
        νe = params.cache.νe[i]
        B  = params.cache.B[i]

        Ωe = e * B / me / νe
        ohmic_heating = me * (1 + Ωe^2) * ue^2 * ne / e * νe
    end

    wall_loss      = ne * params.config.wall_loss_model(U, params, i)
    inelastic_loss = inelastic_losses!(U, params, i)

    return ohmic_heating - wall_loss - inelastic_loss
end