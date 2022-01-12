function apply_reactions!(Q, U, params, i::Int64) #replace Te with Tev
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    reactions, species_range_dict = params.reactions, params.species_range_dict
    _, __, cell_volume = params.z_cell, params.z_edge, params.cell_volume
    dt = params.dt
    index = params.index

    fluid = fluids[1].species.element
    ne = U[index.ne, i]
    neutral_velocity = fluids[1].conservation_laws.u

    for r in reactions
        reactant_index = species_range_dict[r.reactant.symbol][1]
        product_index = species_range_dict[r.product.symbol][1]
        n_reactant = U[reactant_index, i] / fluid.m
        if n_reactant > 1
            k = r.rate_coeff
            @views Q[reactant_index] -= ne * n_reactant * k(U[index.Tev, i]) * fluid.m  #no more *dt/cell_volume, dt taken care of in timemarching scheme, /cell_volume was wrong
                                         #can probably use periodic callback
            @views Q[product_index] += ne * n_reactant * k(U[index.Tev, i]) * fluid.m 
                                        #can probably use periodic callback
            @views Q[product_index + 1] += ne * n_reactant * k(U[index.Tev, i]) * fluid.m *
                                             neutral_velocity #momentum transfer
        end
    end
end

function apply_ion_acceleration!(Q, U, params, i) #make use of calculated potential not electric field input
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    index = params.index
    for j in 1:length(fluids)
        E_d = 0.0
        if i <= length(params.z_cell)
            E_d = -(U[index.ϕ, i] - U[index.ϕ, i - 1]) / (params.z_cell[i] - params.z_cell[i - 1])
        end
        if fluids[j].species.Z > 0
            ni = U[fluid_ranges[j][1], i]
            @views Q[fluid_ranges[j][2]] += e / m(fluids[j]) *
                                            ni *
                                            E_d *
                                            fluids[j].species.Z
        end
    end
end

function source_electron_energy!(Q, U, params, i)
    index = params.index
    uₑ = electron_velocity(U, params, i)
    grad_pe = first_deriv_central_diff(U[index.pe, :], params.z_cell, i)
    grad_ϕ = first_deriv_central_diff(U[index.ϕ, :], params.z_cell, i)
    ν = params.cache.νan[i] + params.cache.νc[i]
    #Hara source term
    #QE = grad_pe*uₑ + mₑ*params.cache.ne[i]*ν*uₑ^2 - S_wall_simple(U[4, :], i) - S_coll(U, params, i) #resistive heating collisions, u has to be total u not just z, azimuthal component dominating
    #Landmark source term
    @views Q[4] = U[index.ne, i]*uₑ*grad_ϕ - S_coll(U, params, i) #U[index.ne, i]*uₑ*grad_ϕ - U[index.ne, i]*S_wall_simple(3/2*U[index.Tev, :], i) - S_coll(U, params, i)
    #=@show Q[4]
    @show i
    @show U[index.ne, i]*uₑ*grad_ϕ
    @show - U[index.ne, i]*S_wall_simple(3/2*U[index.Tev, :], i)
    @show - S_coll(U, params, i)=#
end