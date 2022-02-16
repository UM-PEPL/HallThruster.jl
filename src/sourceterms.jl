function apply_reactions!(Q, U, params, i::Int64) #replace Te with Tev
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    reactions, species_range_dict = params.reactions, params.species_range_dict
    _, __, cell_volume = params.z_cell, params.z_edge, params.cell_volume
    dt = params.dt
    index = params.index

    mi = m(fluids[1])
    ne = U[index.ne, i]
    ϵ = U[index.Tev, i]
    neutral_velocity = fluids[1].conservation_laws.u

    for r in reactions
        reactant_index = species_range_dict[r.reactant.symbol][1]
        product_index = species_range_dict[r.product.symbol][1]
        n_reactant = U[reactant_index, i] / mi
        if n_reactant > 1
            k = r.rate_coeff

            ndot = k(ϵ) * n_reactant * ne
            #can probably use periodic callback
            Q[reactant_index] -= ndot * mi
            Q[product_index] += ndot * mi
            #Q[product_index + 1] += ndot * mi * neutral_velocity #momentum transfer
        end
    end

    # Simplify ionization for now, only single ionization
    #=mi = m(fluids[1])

    k = params.landmark.rate_coeff
    nn = U[1, i] / mi
    ne = U[index.ne, i]
    Te = U[index.Tev, i]

    ndot = k(Te) * nn * ne
    Q[1] -= ndot * mi
    Q[2] += ndot * mi
    Q[3] += ndot * neutral_velocity * mi=#
end

function apply_ion_acceleration!(Q, U, params, i)
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    index = params.index
    for j in 1:length(fluids)
        E_d = -U[index.grad_ϕ, i]
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
    for j in 1:length(fluids)
        if fluids[j].species.Z > 0
            ni = U[fluid_ranges[j][1], i]
            @views Q[fluid_ranges[j][2]] += -e / m(fluids[j]) *
                                            ni * fluids[j].species.Z * U[index.ue, i] / params.cache.μ[i]
        end
    end
end

function source_electron_energy_landmark!(Q, U, params, i)
    index = params.index
    #ν = params.cache.νan[i] + params.cache.νc[i]
    #Hara source term
    #QE = grad_pe*uₑ + mₑ*params.cache.ne[i]*ν*uₑ^2 - S_wall_simple(U[4, :], i) - S_coll(U, params, i) #resistive heating collisions, u has to be total u not just z, azimuthal component dominating
    #Landmark source term
    if params.z_cell[i] <= 0.025
        νε = 1*1e7
    else
        νε = 1*1e7
    end
    UU = 20.0
    W = νε * U[index.Tev, i] * exp(-UU / U[index.Tev, i])
    @views Q[4] =  U[index.ne, i] * (-U[index.ue, i] * -U[index.grad_ϕ, i] - U[1, i]/HallThruster.Xenon.m * params.landmark.loss_coeff(U[index.Tev, i]) - W)
    #@views Q[4] = U[index.ne, i]*uₑ*grad_ϕ - S_coll(U, params, i) - S_wall_simple(U[index.Tev, :], i)
end