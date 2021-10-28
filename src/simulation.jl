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

    return species, fluids, fluid_ranges, species_range_dict
end

function allocate_arrays(sim)
    # Number of variables in the state vector U
    # U = [nn, ni1, ni1ui1..., niN, niNuiN, Te, ne, Φ]
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

	nvariables = size(U, 1)
    ncells = size(U, 2) - 2

	Te_index = nvariables-2
	ne_index = nvariables-1
	ϕ_index = nvariables

    # zero variables and apply boundary conditions
    # TEMPORARY: apply Neumann BC on right edge
    for i in 1:nvariables
	    dU[i, 1] = 0.0
	    dU[i, ncells] = 0.0
        U[i, end] = U[i, end-1]
    end

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

		Δz = z_edge[right] - z_edge[left]

        first_fluid_index = 1
        last_fluid_index = fluid_ranges[end][end]

		for j in first_fluid_index:last_fluid_index
			@views dU[j, i] = (F[j, left] - F[j, right])/Δz + Q[j]
		end
    end
    return nothing
end

left_edge(i) = i-1
right_edge(i) = i

function electron_density(U, fluid_ranges)
    ne = 0.0
    for (i, f) in enumerate(fluid_ranges)
        if i == 1
            continue # neutrals do not contribute to electron density
        end
        charge_state = i-1
        ne += charge_state * U[f[1]]
    end
    return ne
end

function run_simulation(sim)

    species, fluids, fluid_ranges, species_range_dict = configure_simulation(sim)
    z_cell, z_edge = generate_grid(sim.geometry, sim.ncells)

    U, cache = allocate_arrays(sim)

    initial_condition!(U, z_cell, sim, fluid_ranges)

    scheme = sim.scheme

    reactions = load_ionization_reactions(species)

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
    sol = solve(prob, SSPRK33(), dt = sim.dt, saveat = [])
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

    Te_index = nvariables - 2
    ne_index = nvariables - 1
    ϕ_index  = nvariables

    for (i, z) in enumerate(z_cell)
        U[nn_index, i] = nn

        ni = sim.initial_ni(z)
        Te = sim.initial_Te(z)
        ϕ = sim.initial_ϕ(z)

        # ions initialized with equal densities, same velocity as neutrals
        for j in fluid_ranges[2:end]
            n_index = j[1]
            nu_index = j[2]
            U[n_index, i] = ni
            U[nu_index, i] = ni * un
        end
        U[Te_index, i] = Te
        @views U[ne_index, i] = electron_density(U[:, i], fluid_ranges)
        U[ϕ_index, i]  = ϕ
    end

    return U
end
