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
            ni = U[fluid_ranges[j][1], i]
            @views Q[fluid_ranges[j][2]] += e/m(fluids[j])*ni*E_d*fluids[j].species.Z
        end
    end
end