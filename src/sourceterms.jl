function apply_reactions!(Q, U, params, i::Int64)
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
	reactions, species_range_dict = params.reactions, params.species_range_dict
    _, __, cell_volume = params.z_cell, params.z_edge, params.cell_volume

    fluid = fluids[1].species.element
	ne = @views HallThruster.electron_density(U[:, i], fluid_ranges)/fluid.m
    Te = 15.0 #in eV
    dt = 1e-6
	for r in reactions
		reactant_index = species_range_dict[r.reactant][1]
		product_index = species_range_dict[r.product][1]
		n_reactant = U[reactant_index, i]/fluid.m
		if n_reactant > 1
			n_product = U[product_index, i]/fluid.m
			k = r.rate_coeff
			@views Q[reactant_index] -= ne * n_reactant * k(Te) * dt*fluid.m/cell_volume #can probably use periodic callback
			println("$(Q[reactant_index])")
			@views Q[product_index]  += ne * n_product  * k(Te) * dt*fluid.m/cell_volume #can probably use periodic callback
			println("$(Q[product_index])")
		end
	end
end

function apply_ion_acceleration!(Q, U, params, i)
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    for j in 1:length(fluids)
        if fluids[j].species.Z > 0
            @views Q[fluid_ranges[j][2]] = e/fluids[j].species.element.m*U[fluid_ranges[j][1], i]*params.E_d[i]*fluids[j].species.Z
        end
    end
end