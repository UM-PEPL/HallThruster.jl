# Analytic radiative decay transfers heavy-species mass and momentum.
# It does not affect electrons or constrain the timestep.

function apply_deexcitation_reactions!(fluid_arr, params)
    (;
        deexcitation_reactions,
        deexcitation_reactant_indices,
        deexcitation_product_indices,
        landmark,
    ) = params

    # Use the timestep applied by SSPRK, not the next CFL estimate.
    dt = params.dt[]

    rxns = zip(
        deexcitation_reactions, deexcitation_reactant_indices, deexcitation_product_indices,
    )

    for (rxn, reactant_index, product_index) in rxns
        apply_deexcitation!(fluid_arr, reactant_index, product_index, rxn.rates, dt, landmark)
    end
    return Inf
end

function apply_deexcitation!(fluids, reactant_index, product_index, rates, dt, landmark)
    reactant = fluids[reactant_index]
    inv_m = 1 / reactant.species.element.m
    ncells = length(reactant.density)
    total_rate = sum(rates)
    (total_rate > 0 && dt > 0) || return Inf

    # Preserve branch ratios while integrating total decay over the step.
    frac = -expm1(-total_rate * dt)
    eff_scale = frac / (total_rate * dt)

    # Reactant loss
    effective_rate = eff_scale * total_rate
    @inbounds @simd for i in 2:(ncells - 1)
        reactant.dens_ddt[i] -= effective_rate * reactant.density[i]

        if !landmark && reactant.type != _ContinuityOnly
            reactant.mom_ddt[i] -= effective_rate * reactant.momentum[i]
        end
    end

    # Product gain
    @inbounds for (prod_ind, rate) in zip(product_index, rates)
        product = fluids[prod_ind]
        mass_ratio = product.species.element.m * inv_m
        branch_rate = eff_scale * rate

        @simd for i in 2:(ncells - 1)
            product.dens_ddt[i] += mass_ratio * branch_rate * reactant.density[i]

            if !landmark && reactant.type != _ContinuityOnly && product.type != _ContinuityOnly
                product.mom_ddt[i] += mass_ratio * branch_rate * reactant.momentum[i]
            end
        end
    end

    return Inf
end
