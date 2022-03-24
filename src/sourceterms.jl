function apply_reactions!(Q, U, params, i::Int64) #replace Te with Tev
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    reactions, species_range_dict = params.reactions, params.species_range_dict
    (;index, fluids, fluid_ranges, reactions, species, z_cell, z_edge, cell_volume, dt, index) = params
    dt = params.dt
    index = params.index

    mi = m(fluids[1])
    ne = params.cache.ne[i]
    ϵ = U[index.nϵ] / ne

    for r in reactions
        reactant_index = species_range_dict[r.reactant.symbol][1]
        product_index = species_range_dict[r.product.symbol][1]
        n_reactant = U[reactant_index, i] / mi
        if n_reactant > 1
            k = r.rate_coeff

            ndot = k(ϵ) * n_reactant * ne
            Q[reactant_index] -= ndot * mi
            Q[product_index] += ndot * mi
        end
    end
end

function apply_ion_acceleration!(Q, U, params, i)
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    index = params.index
    ∇ϕ = params.cache.∇ϕ
    for j in 1:length(fluids)
        E_d = -∇ϕ[i]
        if fluids[j].species.Z > 0
            ni = U[fluid_ranges[j][1], i]
            @views Q[fluid_ranges[j][2]] += e / m(fluids[j]) *
                                            ni *
                                            E_d *
                                            fluids[j].species.Z
        end
    end
end

function apply_ion_acceleration_coupled!(Q, U, params, i)
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    index = params.index
    (;ue, μ) = params.cache
    for j in 1:length(fluids)
        if fluids[j].species.Z > 0
            ni = U[fluid_ranges[j][1], i]
            @views Q[fluid_ranges[j][2]] += -e / m(fluids[j]) *
                                            ni * fluids[j].species.Z * ue[i] / μ[i]
        end
    end
end


function source_electron_energy!(Q, U, params, i)
    Q[params.index.nϵ] = source_electron_energy(U, params, i)
end

function source_electron_energy(U, params, i)
    index = params.index

    mi = params.propellant.m
    ne = params.cache.ne[i]
    ϵ = U[index.nϵ, i] / ne
    ue = params.cache.ue[i]
    ∇ϕ = params.cache.∇ϕ[i]

    nn = U[index.ρn, i] / mi
    K = params.collisional_loss_coeff(ϵ)
    W = params.config.wall_loss_model(U, params, i)

    ohmic_heating      = ne * ue * ∇ϕ
    wall_losses        = ne * W
    collisional_losses = ne * nn * K

    return ohmic_heating - wall_losses - collisional_losses
end