function update!(du, u, params, t)
	fluids, fluid_ranges = params.fluids, params.fluid_ranges
	reactions, species_range_dict = params.reactions, params.species_range_dict

	F, UL, UR, Q = params.cache

	z_cell, z_edge = params.z_cell, params.z_edge
	scheme = params.scheme

	nvariables, ncells = size(U, 2)

	ne_index = nvariables-2
	Te_index = nvariables-1
	ϕ_index = nvariables

	nedges = size(F, 2)

	dU[:, 1] .= 0.0
	dU[:, ncells] .= 0.0

	reconstruct!(UL, UR, U, scheme)

	# Compute heavy species fluxes
    for i in 1:nedges
		for (fluid, fluid_range) in zip(fluids, fluid_ranges)
			@views F[fluid_range, i] .= scheme.flux_function(
				UL[fluid_range, i],
				UR[fluid_range, i],
				fluid
			)
		end
	end

	# Compute heavy species source terms
	for i in 2:ncells-1

		Q .= 0.0

		# Compute heavy species source term due to electric field
		for (fluid, fluid_range) in zip(fluids, fluid_ranges)

			if fluid.species.Z == 0
				continue # Neutrals not affected by electric field
			end
			density_index = fluid_range[1]
			momentum_index = fluid_range[2]
			ΔΦ = U[ϕ_index, i+1] - U[ϕ_index, i-1]
			Δz = z_cell[i+1] -  z_cell[i-1]
			E = -ΔΦ / Δz
			q = e * fluid.species.Z
			n = U[density_index]
			Q[momentum_index] += q * n * E / m(fluid)
		end

		# Compute electron density in cell
		ne = @views electron_density(U[:, i], fluid_ranges)
		U[ne_index, i] = ne

		Te = U[Te_index, i]

		# Compute heavy species source term due to ionization
		for r in reactions
			reactant_index = species_range_dict[r.reactant][1]
			product_index = species_range_dict[r.product][1]
			n_reactant = U[reactant_index, i]
			n_product = U[product_index, i]
			k = r.rate_coeff
			Q[reactant_index] -= ne * n_reactant * k(Te)
			Q[product_index]  += ne * n_product  * k(Te)
		end

		# Compute dU/dt
		left = left_edge(i)
		right = right_edge(i)

		Δx = x_edge[right] - x_edge[left]

		for j in (fluid_ranges[1][1]):(fluid_ranges[end][end])
			@views dU[j, i] = (F[j, left] - F[j, right])/Δx + Q[j, i]
		end
    end
end

left_edge(i) = i-1
right_edge(i) = i

electron_density(u, fluid_ranges) = sum(u[r[1]] for r in fluid_ranges)


