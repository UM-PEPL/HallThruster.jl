function apply_reactions!(Q, U, params, Tev, i::Int64) #replace Te with Tev
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
	reactions, species_range_dict = params.reactions, params.species_range_dict
    _, __, cell_volume = params.z_cell, params.z_edge, params.cell_volume
    dt = params.dt

    fluid = fluids[1].species.element
	ne = @views electron_density(U[:, i], fluid_ranges)/fluid.m
    neutral_velocity = fluids[1].conservation_laws.u
    
	for r in reactions
		reactant_index = species_range_dict[r.reactant.symbol][1]
		product_index = species_range_dict[r.product.symbol][1]
		n_reactant = U[reactant_index, i]/fluid.m
		if n_reactant > 1
			k = r.rate_coeff
			@views Q[reactant_index] -= ne * n_reactant * k(Tev[i]) * dt*fluid.m/cell_volume #can probably use periodic callback
			@views Q[product_index]  += ne * n_reactant  * k(Tev[i]) * dt*fluid.m/cell_volume #can probably use periodic callback
            @views Q[product_index+1] += ne * n_reactant  * k(Tev[i]) * dt*fluid.m/cell_volume * neutral_velocity #momentum transfer
		end
	end
end

function apply_ion_acceleration!(Q, U, params, ϕ, i) #make use of calculated potential not electric field input
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    for j in 1:length(fluids)
        E_d = 0.0
        if i < 101
            E_d = -(ϕ[i]-ϕ[i-1])/(params.z_cell[i] - params.z_cell[i-1])
        end
        if fluids[j].species.Z > 0
            @views Q[fluid_ranges[j][2]] += e/fluids[j].species.element.m*U[fluid_ranges[j][1], i]*E_d*fluids[j].species.Z
        end
    end
end

function set_up_potential_equation!(U, A, b, Tev, params)
    #directly discretising the equation, conserves properties such as negative semidefinite etc... 
    #add functionality for nonuniform cell size
    z_cell, z_edge, _ = params.z_cell, params.z_edge, params.cell_volume
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    fluid = fluids[1].species.element
    N = length(z_cell) - 2 #remove boundary ghost cells

    #have to think about how to do this, either with a vector or a function that will be called multiple times using Tev and electron density
    p = params.cache.pe
    B = params.cache.B

    for i in 1:N+2
        @views p[i] = electron_pressure(electron_density(U[:, i], fluid_ranges)/fluid.m, Tev[i])
    end

    #using Dirichlet boundary conditions
    #make smth like a BC function for this

    ϕ_L = 400
    ϕ_R = 0

    Δz = z_edge[3] - z_edge[2]
    function source_term_potential() #make this to be able to verify implementation with 0 source term
    end

    #left boundary
    @views ne⁻ = electron_density(U[:, 1], fluid_ranges)/fluid.m #on boundary, ghost cell, perfect fit for i-1/2 on stencil 
    @views ne⁺ = (electron_density(U[:, 2], fluid_ranges)/fluid.m + electron_density(U[:, 3], fluid_ranges)/fluid.m)/2
    nn⁻ = U[1, 1]/fluid.m
    nn⁺ = (U[1, 2]/fluid.m + U[1, 3]/fluid.m)/2
    μ⁻ = cf_electron_transport(get_v_an(), get_v_c(Tev[1], ne⁻, nn⁻, fluid.m), B[1])#B_field(B_max, z_cell[1], L_ch))
    μ⁺ = cf_electron_transport(get_v_an(), get_v_c((Tev[2] + Tev[3])/2, ne⁺, nn⁺, fluid.m), 0.5 * (B[2] + B[3]))#B_field(B_max, (z_cell[3]+z_cell[2])/2, L_ch))
    #A[1, 1] = -1.5*(ne⁻*μ⁻ + ne⁺*μ⁺)/(Δz)^2 
    A[1, 1] = -1*(ne⁻*μ⁻ + ne⁺*μ⁺)/(Δz)^2 #no interpolation on boundary
    A[1, 2] = ne⁺*μ⁺/(Δz)^2
    #b[1] = - ϕ_L*(ne⁻*μ⁻ + ne⁺*μ⁺)/((Δz)^2) #- 1.5*(μ⁺ + μ⁻)*p[2]/(Δz)^2 + μ⁺*p[3]/(Δz)^2 + (μ⁺ + μ⁻)*p[1]/(Δz)^2 +
    #- U[3, 1]/(Δz) + (U[3, 3] + U[3, 2])/(2*Δz)
    b[1] = - ϕ_L*(ne⁻*μ⁻)/((Δz)^2) - 1*(μ⁺ + μ⁻)*p[2]/(Δz)^2 + μ⁺*p[3]/(Δz)^2 + (μ⁻)*p[1]/(Δz)^2
    - U[3, 1]/(Δz) + (U[3, 3] + U[3, 2])/(2*Δz) #no interpolation on boundary

    #right boundary
    @views ne⁻ = (electron_density(U[:, N+1], fluid_ranges)/fluid.m + electron_density(U[:, N], fluid_ranges)/fluid.m)/2
    @views ne⁺ =  electron_density(U[:, N+2], fluid_ranges)/fluid.m
    nn⁻ = (U[1, N+1]/fluid.m + U[1, N]/fluid.m)/2
    nn⁺ = U[1, N+2]/fluid.m
    μ⁻ = cf_electron_transport(get_v_an(), get_v_c((Tev[N+1] + Tev[N])/2, ne⁻, nn⁻, fluid.m), 0.5 * (B[N] + B[N+1]))#B_field(B_max, (z_cell[N]+z_cell[N+1])/2, L_ch))
    μ⁺ = cf_electron_transport(get_v_an(), get_v_c(Tev[N+2], ne⁺, nn⁺, fluid.m), B[N+2])#B_field(B_max, z_cell[N+2], L_ch))
    #A[N, N] = -1.5*(ne⁻*μ⁻ + ne⁺*μ⁺)/(Δz)^2
    A[N, N] = -1*(ne⁻*μ⁻ + ne⁺*μ⁺)/(Δz)^2 #no interpolation on boundary
    A[N, N-1] = ne⁻*μ⁻/(Δz)^2
    #b[N] = - ϕ_R*(ne⁻*μ⁻ + ne⁺*μ⁺)/((Δz)^2) #+ μ⁻*p[N]/(Δz)^2 - 1.5*(μ⁺ + μ⁻)*p[N+1]/(Δz)^2 + (μ⁺ + μ⁻)*p[N+2]/(Δz)^2
    #- (U[3, N] + U[3, N+1])/(2*Δz) + U[3, N+2]/(Δz)
    b[N] = - ϕ_R*(ne⁺*μ⁺)/((Δz)^2) + μ⁻*p[N]/(Δz)^2 - 1*(μ⁺ + μ⁻)*p[N+1]/(Δz)^2 + (μ⁺)*p[N+2]/(Δz)^2
    - (U[3, N] + U[3, N+1])/(2*Δz) + U[3, N+2]/(Δz) #no interpolation on boundary

    for i in 2:N-1
        i_f = i+1 #fluid index, due to ghost cell on boundary
        @views ne⁻ = (electron_density(U[:, i_f], fluid_ranges)/fluid.m + electron_density(U[:, i_f-1], fluid_ranges)/fluid.m)/2
        @views ne⁺ = (electron_density(U[:, i_f+1], fluid_ranges)/fluid.m + electron_density(U[:, i_f], fluid_ranges)/fluid.m)/2
        nn⁻ = (U[1, i_f]/fluid.m + U[1, i_f-1]/fluid.m)/2
        nn⁺ = (U[1, i_f]/fluid.m + U[1, i_f+1]/fluid.m)/2
        μ⁻ = cf_electron_transport(get_v_an(), get_v_c((Tev[i_f] + Tev[i_f-1])/2, ne⁺, nn⁺, fluid.m), 0.5 * (B[i_f-1] + B[i_f]))#B_field(B_max, (z_cell[i_f-1]+z_cell[i_f])/2, L_ch))
        μ⁺ = cf_electron_transport(get_v_an(), get_v_c((Tev[i_f] + Tev[i_f+1])/2, ne⁺, nn⁺, fluid.m), 0.5 * (B[i_f] + B[i_f+1]))#B_field(B_max, (z_cell[i_f+1]+z_cell[i_f])/2, L_ch))
		Δz = z_edge[i-1] - z_edge[i]
        #=Δz¹ = z_edge[i-2] - z_edge[i-1]
        Δz² = z_edge[i-1] - z_edge[i]=#
        A[i, i-1] = ne⁻*μ⁻/(Δz)^2
        A[i, i] = -(ne⁻*μ⁻ + ne⁺*μ⁺)/(Δz)^2
        A[i, i+1] = ne⁺*μ⁺/(Δz)^2
        b[i] = μ⁻*p[i_f-1]/(Δz)^2 - (μ⁺ + μ⁻)*p[i_f]/(Δz)^2 + μ⁺*p[i_f+1]/(Δz)^2
        + U[3, i_f+1]/(2*Δz) - U[3, i_f-1]/(2*Δz)
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