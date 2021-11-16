function apply_reactions!(Q, U, params, i::Int64)
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
	reactions, species_range_dict = params.reactions, params.species_range_dict
    _, __, cell_volume = params.z_cell, params.z_edge, params.cell_volume
    dt = params.dt

    fluid = fluids[1].species.element
	ne = @views HallThruster.electron_density(U[:, i], fluid_ranges)/fluid.m
    Te = 2.0 #in eV
    neutral_velocity = fluids[1].conservation_laws.u
    
	for r in reactions
		reactant_index = species_range_dict[r.reactant][1]
		product_index = species_range_dict[r.product][1]
		n_reactant = U[reactant_index, i]/fluid.m
		if n_reactant > 1
			k = r.rate_coeff
			@views Q[reactant_index] -= ne * n_reactant * k(Te) * dt*fluid.m/cell_volume #can probably use periodic callback
			@views Q[product_index]  += ne * n_reactant  * k(Te) * dt*fluid.m/cell_volume #can probably use periodic callback
            @views Q[product_index+1] += ne * n_reactant  * k(Te) * dt*fluid.m/cell_volume * neutral_velocity #momentum transfer
		end
	end
end

function apply_ion_acceleration!(Q, U, params, i)
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    for j in 1:length(fluids)
        if fluids[j].species.Z > 0
            @views Q[fluid_ranges[j][2]] += e/fluids[j].species.element.m*U[fluid_ranges[j][1], i]*params.E_d[i]*fluids[j].species.Z
        end
    end
end