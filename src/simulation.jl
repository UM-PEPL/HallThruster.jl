get_species(sim) = [Species(sim.propellant, i) for i in 0:sim.ncharge]

function configure_simulation(sim)
    species = get_species(sim)
    fluids = [
        Fluid(species[1], ContinuityOnly(sim.neutral_velocity, sim.neutral_temperature));
        [Fluid(species[i], IsothermalEuler(sim.ion_temperature)) for i in 2:sim.ncharge+1]
    ]
    fluid_ranges = ranges(fluids)
    species_range_dict = Dict(
        fluid.species => fluid_range for (fluid, fluid_range) in zip(fluids, fluid_ranges)
    )

    return fluids, fluid_ranges, species_range_dict
end

function allocate_arrays(sim)
    # Number of variables in the state vector U
    # U = [nn, ni1, ni1, ni1ui1..., niN, niNuiN, Te, ne, Φ]
    nvariables = 1 + 2 * sim.ncharge + 3

    ncells = sim.ncells
    nedges = sim.ncells + 1

    U = zeros(nvariables, ncells+2) # need to allocate room for ghost cells
    F = zeros(nvariables, nedges)
    UL = zeros(nvariables, nedges)
    UR = zeros(nvariables, nedges)
    Q = zeros(nvariables)

    cache = (F, UL, UR, Q)
    return U, cache
end

function update!(dU, U, params, t)
	fluids, fluid_ranges = params.fluids, params.fluid_ranges
	reactions, species_range_dict = params.reactions, params.species_range_dict

	F, UL, UR, Q = params.cache

	z_cell, z_edge = params.z_cell, params.z_edge
	scheme = params.scheme

	nvariables, ncells = size(U, 2)

	ne_index = nvariables-2
	Te_index = nvariables-1
	ϕ_index = nvariables

	dU[:, 1] .= 0.0
	dU[:, ncells] .= 0.0

	# TEMPORARY: apply Neumann BC on right edge
	U[:, end] .= U[:, end-1]

	reconstruct!(UL, UR, U, scheme)
	compute_fluxes!(F, UL, UR, fluids, fluid_ranges, scheme)

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

function run_simulation(sim)

    fluids, fluid_ranges, species_range_dict = configure_simulation(sim)
    z_cell, z_edge = generate_grid(sim.geometry, sim.ncells)

    U, cache = allocate_arrays(sim)

    initial_condition!(U, z_cell, sim, fluid_ranges)

    scheme = sim.scheme
    reactions = sim.reactions

    params = (;
        cache,
        fluids,
        fluid_ranges,
        species_range_dict,
        z_cell,
        z_edge,
        reactions,
        scheme
    )

    prob = ODEProblem{true}(update!, U, sim.tspan, params)
    sol = solve(prob, saveat = [])
    return sol
end

function inlet_neutral_density(sim)
    un = sim.neutral_velocity
    A = channel_area(sim.geometry)
    m_atom = sim.propellant.m
    nn = sim.inlet_mdot / un / A / m_atom
    return nn
end

function initial_condition!(U, z_cell, sim, fluid_ranges)
    nvariables = size(U, 1)
    nn = inlet_neutral_density(sim)
    un = sim.neutral_velocity

    nn_index = 1
    ni_index = 2
    ni_ui_index = 3

    Te_index = nvariables - 2
    ne_index = nvariables - 1
    ϕ_index  = nvariables

    for (i, z) in enumerate(z_cell)
        U[nn_index, i] = nn

        ne = sim.initial_ne(z)
        Te = sim.initial_Te(z)
        ϕ = sim.initial_ϕ(z)

        # singly-charged ions initialized with same density as electrons,
        # and same velocity as neutrals
        U[ni_index, i] = ne
        U[ni_ui_index, i] = ne * un

        for j in fluid_ranges[3:end]
            # Initialize all other charge states to zero
            U[j, i] .= 0.0
        end
        U[Te_index, i] = Te
        U[ne_index, i] = ne
        U[ϕ_index, i]  = ϕ
    end

    return U
end
