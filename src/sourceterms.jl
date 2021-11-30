function apply_reactions!(Q, U, params, Tev, i::Int64) #replace Te with Tev
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    reactions, species_range_dict = params.reactions, params.species_range_dict
    _, __, cell_volume = params.z_cell, params.z_edge, params.cell_volume
    dt = params.dt

    fluid = fluids[1].species.element
    ne = @views electron_density(U[:, i], fluid_ranges) / fluid.m
    neutral_velocity = fluids[1].conservation_laws.u

    for r in reactions
        reactant_index = species_range_dict[r.reactant.symbol][1]
        product_index = species_range_dict[r.product.symbol][1]
        n_reactant = U[reactant_index, i] / fluid.m
        if n_reactant > 1
            k = r.rate_coeff
            @views Q[reactant_index] -= ne * n_reactant * k(Tev[i]) * dt * fluid.m /
                                        cell_volume #can probably use periodic callback
            @views Q[product_index] += ne * n_reactant * k(Tev[i]) * dt * fluid.m /
                                       cell_volume #can probably use periodic callback
            @views Q[product_index + 1] += ne * n_reactant * k(Tev[i]) * dt * fluid.m /
                                           cell_volume * neutral_velocity #momentum transfer
        end
    end
end

function apply_ion_acceleration!(Q, U, params, ϕ, i) #make use of calculated potential not electric field input
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    for j in 1:length(fluids)
        E_d = 0.0
        if i < length(params.z_cell) - 1
            E_d = -(ϕ[i] - ϕ[i - 1]) / (params.z_cell[i] - params.z_cell[i - 1])
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

#=
function set_up_potential_equation_staggered!(U, A, b, Tev, params)
    #checkerboard
    #add functionality for nonuniform cell size
    z_cell, z_edge, _ = params.z_cell, params.z_edge, params.cell_volume
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    fluid = fluids[1].species.element
    N = length(z_cell-2) #remove boundary ghost cells

    #have to think about how to do this, either with a vector or a function that will be called multiple times using Tev and electron density
    p = zeros(N)
    p .= electron_pressure(10e19, Tev)

    #cross field electron transport constant for now
    #would change with B, neutral density, and depending on model for classical and anomalous transport more variables
    cf_transport⁻ = cf_electron_transport(get_v_an(), get_v_c(), B_field())
    μ⁺ = cf_electron_transport(get_v_an(), get_v_c(), B_field())

    #need to define electron density and cf transport for each cell on boundaries. 
    #using Dirichlet boundary conditions
    #make smth like a BC function for this
    #replace p[i] on RHS with boundary state, can do the same with cf transport. Can use boundary faces (ghost cells)
    #which are included in the fluid solve

    ϕ_L = 450
    ϕ_R = 10
    #left boundary
    @views A[1, 1] = -2.5*(ne⁻*cf_transport⁻ + ne⁺*μ⁺)/(2*Δz)^2
    @views A[1, 3] = ne⁺*μ⁺/(2*Δz)^2
    @views b[1] = -2.5*(μ⁺ + cf_transport⁻)*p[i]/(2*Δz)^2 + μ⁺*p[i+2]/(2*Δz)^2 + 2*(μ⁺ + cf_transport⁻)*p[i]/(2*Δz)^2 +
    + U[3, i-1]/(2*Δz) + U[3, i+1]/(2*Δz) - 4*ϕ_L/((2*Δz)^2)

    @views A[2, 2] = -7/6*(ne⁻*cf_transport⁻ + ne⁺*μ⁺)/(2*Δz)^2
    @views A[2, 4] = ne⁺*μ⁺/(2*Δz)^2
    @views b[2] = -7/6*(μ⁺ + cf_transport⁻)*p[i]/(2*Δz)^2 + μ⁺*p[i+2]/(2*Δz)^2 + 2/3*(μ⁺ + cf_transport⁻)*p[i]/(2*Δz)^2
    + U[3, i-1]/(2*Δz) + U[3, i+1]/(2*Δz) - 4/3*ϕ_L/((2*Δz)^2)

    #right boundary
    @views A[N-1, N-1] = -7/6*(ne⁻*cf_transport⁻ + ne⁺*μ⁺)/(2*Δz)^2
    @views A[N-1, N-3] = ne⁻*cf_transport⁻/(2*Δz)^2
    @views b[N-1] = cf_transport⁻*p[i-2]/(2*Δz)^2 - 7/6*(μ⁺ + cf_transport⁻)*p[i]/(2*Δz)^2 + 2/3*(μ⁺ + cf_transport⁻)*p[i]/(2*Δz)^2
    + U[3, i-1]/(2*Δz) + U[3, i+1]/(2*Δz) - 4/3*ϕ_R/((2*Δz)^2)

    @views A[N, N] = -2.5*(ne⁻*cf_transport⁻ + ne⁺*μ⁺)/(2*Δz)^2
    @views A[N, N-2] = ne⁻*cf_transport⁻/(2*Δz)^2
    @views b[N] = cf_transport⁻*p[i-2]/(2*Δz)^2 - 2.5*(μ⁺ + cf_transport⁻)*p[i]/(2*Δz)^2 + 2*(μ⁺ + cf_transport⁻)*p[i]/(2*Δz)^2
    + U[3, i-1]/(2*Δz) + U[3, i+1]/(2*Δz) - 4*ϕ_R/((2*Δz)^2)

    for i in 3:length(z_cell)-2
        #cf_transport⁻ = cf_electron_transport(get_v_an(), get_v_c(), B_field())
        #μ⁺ = cf_electron_transport(get_v_an(), get_v_c(), B_field())
        ne⁻ = electron_density(U[:, i-1], fluid_ranges)/fluid.m
        ne⁺ = electron_density(U[:, i+1], fluid_ranges)/fluid.m
    
		Δz = z_edge[i-3] - z_edge[i-2]
        #=Δz¹ = z_edge[i-2] - z_edge[i-1]
        Δz² = z_edge[i-1] - z_edge[i]
        Δz³ = z_edge[i] - z_edge[i+1]
        Δz⁴ = z_edge[i+1] - z_edge[i+2]=#
        @views A[i, i-2] = ne⁻*cf_transport⁻/(2*Δz)^2
        @views A[i, i] = -(ne⁻*cf_transport⁻ + ne⁺*μ⁺)/(2*Δz)^2
        @views A[i, i+2] = ne⁺*μ⁺/(2*Δz)^2
        #assuming there is a vector for electron pressure, will have one for T, need to calculate p from there
        #need to somewhere implement electron pressure from electron temperature, actually want it in eV due to ionization unit
        @views b[i] = cf_transport⁻*p[i-2]/(2*Δz)^2 - (μ⁺ + cf_transport⁻)*p[i]/(2*Δz)^2 + μ⁺*p[i+2]/(2*Δz)^2
        + U[3, i-1]/(2*Δz) + U[3, i+1]/(2*Δz)
    end 

end 

=#