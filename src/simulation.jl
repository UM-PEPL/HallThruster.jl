struct HyperbolicScheme{F, L}
    flux_function::F  # in-place flux function
    limiter::L # limiter
    reconstruct::Bool
end

Base.@kwdef mutable struct MultiFluidSimulation{IC, B1, B2, S, F, L} #could add callback, or autoselect callback when in MMS mode
    grid::Grid1D
    fluids::Vector{Fluid}     # An array of user-defined fluids. 
                              # This will give us the capacity to more easily do shock tubes (and other problems)
                              # without Hall thruster baggage
    initial_condition::Vector{IC}
    boundary_conditions::Tuple{B1, B2}   # Tuple of left and right boundary conditions, subject to the approval of PR #10 
    end_time::Float64    # How long to simulate
    scheme:: HyperbolicScheme{F, L} # Flux, Limiter
    source_term!::S  # Source term function. This can include reactons, electric field, and MMS terms
    saveat::Vector{Float64} #when to save
end

get_species(sim) = [Species(sim.propellant, i) for i in 0:sim.ncharge]

function configure_simulation(sim)
    species = get_species(sim)
    fluids = [
        Fluid(species[1], ContinuityOnly(sim.neutral_velocity, sim.neutral_temperature));
        #[Fluid(species[i], IsothermalEuler(sim.ion_temperature)) for i in 2:sim.ncharge+1]
        [Fluid(species[i], EulerEquations()) for i in 2:sim.ncharge+1]
    ]
    fluid_ranges = ranges(fluids)
    species_range_dict = Dict(
        fluid.species => fluid_range for (fluid, fluid_range) in zip(fluids, fluid_ranges)
    )

    return species, fluids, fluid_ranges, species_range_dict
end

function configure_simulation_v2(sim)
    fluids = sim.fluids
    species = [fluids[i].species for i in 1:length(fluids)]
    fluid_ranges = ranges(fluids)
    species_range_dict = Dict(
        fluid.species => fluid_range for (fluid, fluid_range) in zip(fluids, fluid_ranges)
    )

    return species, fluids, fluid_ranges, species_range_dict
end

function allocate_arrays(sim)
    # Number of variables in the state vector U
    # U = [nn, ni1, ni1ui1..., niN, niNuiN, Te, ne, Φ]
    #nvariables = 1 + 2 * sim.ncharge + 3
    nvariables = 1 + 3 * sim.ncharge + 3

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

function allocate_arrays_v2(sim) #rewrite allocate arrays as function of set of equations, either 1, 2 or 3
    # Number of variables in the state vector U
    nvariables = 0
    for i in 1:length(sim.fluids)
        if sim.fluids[i].conservation_laws.type == :ContinuityOnly
            nvariables += 1
        elseif sim.fluids[i].conservation_laws.type ==  :IsothermalEuler
            nvariables += 2
        elseif sim.fluids[i].conservation_laws.type == :EulerEquations
            nvariables += 3
        end
    end
    
    ncells = sim.grid.ncells
    nedges = sim.grid.ncells + 1

    U = zeros(nvariables, ncells+2) # need to allocate room for ghost cells
    F = zeros(nvariables, nedges)
    UL = zeros(nvariables, nedges)
    UR = zeros(nvariables, nedges)
    Q = zeros(nvariables)

    cache = (F, UL, UR, Q)
    return U, cache
end

function update_v2!(dU, U, params, t)
	fluids, fluid_ranges = params.fluids, params.fluid_ranges
	reactions, species_range_dict = params.reactions, params.species_range_dict

	F, UL, UR, Q = params.cache

	z_cell, z_edge = params.z_cell, params.z_edge
	scheme = params.scheme

	nvariables = size(U, 1)
    ncells = size(U, 2) - 2

    apply_bc!(U, params.BCs[1], :left)
    apply_bc!(U, params.BCs[2], :right)

    reconstruct!(UL, UR, U, scheme)
	compute_fluxes!(F, UL, UR, fluids, fluid_ranges, scheme)

	# Compute heavy species source terms
	for i in 2:ncells+1 #+1 since ncells takes the amount of cells, but there are 2 more face values
		Q .= 0.0
        #add possible source term input here

		# Compute dU/dt
		left = left_edge(i)
		right = right_edge(i)

		Δz = z_edge[right] - z_edge[left]

        first_fluid_index = 1
        last_fluid_index = fluid_ranges[end][end]

        @views @. dU[:, i] = (F[:, left] - F[:, right])/Δz + Q
    end
    return nothing
end


function update!(dU, U, params, t)
	fluids, fluid_ranges = params.fluids, params.fluid_ranges
	reactions, species_range_dict = params.reactions, params.species_range_dict

	F, UL, UR, Q = params.cache
    MMS, mms! = params.mms

	z_cell, z_edge = params.z_cell, params.z_edge
	scheme = params.scheme

	nvariables = size(U, 1)
    ncells = size(U, 2) - 2

	#=Te_index = nvariables-2
	ne_index = nvariables-1
	ϕ_index = nvariables=#

    apply_bc!(U, params.BCs[1], :left)
    apply_bc!(U, params.BCs[2], :right)

    reconstruct!(UL, UR, U, scheme)
	compute_fluxes!(F, UL, UR, fluids, fluid_ranges, scheme)

	# Compute heavy species source terms
	for i in 2:ncells+1 #+1 since ncells takes the amount of cells, but there are 2 more face values
		Q .= 0.0
        if MMS
            mms!(Q, [z_cell[i]])#mms!(Q, [z_cell[i]])
        end

		# Compute heavy species source term due to electric field
		#=or (fluid, fluid_range) in zip(fluids, fluid_ranges)
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
		end=#

		# Compute electron density in cell
		#ne = @views electron_density(U[:, i], fluid_ranges)
		#U[ne_index, i] = ne

		#Te = U[Te_index, i]

		# Compute heavy species source term due to ionization
		#=for r in reactions
			reactant_index = species_range_dict[r.reactant][1]
			product_index = species_range_dict[r.product][1]
			n_reactant = U[reactant_index, i]
			n_product = U[product_index, i]
			k = r.rate_coeff
			Q[reactant_index] -= ne * n_reactant * k(Te)
			Q[product_index]  += ne * n_product  * k(Te)
		end=#

		# Compute dU/dt
		left = left_edge(i)
		right = right_edge(i)

		Δz = z_edge[right] - z_edge[left]

        first_fluid_index = 1
        last_fluid_index = fluid_ranges[end][end]

        @views @. dU[:, i] = (F[:, left] - F[:, right])/Δz + Q
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
    grid = generate_grid(sim.geometry, sim.ncells)

    U, cache = allocate_arrays(sim)
    mms = (sim.MMS, sim.mms!)

    initial_condition!(U, grid.cell_centers, sim, fluid_ranges)

    scheme = sim.scheme

    reactions = load_ionization_reactions(species)
    BCs = sim.BCs

    params = (;
        cache,
        fluids,
        fluid_ranges,
        species_range_dict,
        z_cell = grid.cell_centers,
        z_edge = grid.edges,
        reactions,
        scheme,
        mms,
        BCs
    )

    prob = ODEProblem{true}(update!, U, sim.tspan, params)
    sol = solve(prob, Tsit5(), saveat = sim.saveat)
    return sol
end

function run_simulation_v2(sim)

    species, fluids, fluid_ranges, species_range_dict = configure_simulation_v2(sim)
    grid = sim.grid

    U, cache = allocate_arrays_v2(sim)

    initial_condition_v2!(U, grid.cell_centers, sim.initial_condition, fluid_ranges)

    scheme = sim.scheme
    tspan = (0., sim.end_time)

    reactions = load_ionization_reactions(species)
    BCs = sim.boundary_conditions

    params = (;
        cache,
        fluids,
        fluid_ranges,
        species_range_dict,
        z_cell = grid.cell_centers,
        z_edge = grid.edges,
        reactions,
        scheme,
        BCs
    )

    prob = ODEProblem{true}(update_v2!, U, tspan, params)
    sol = solve(prob, Tsit5(), saveat = sim.saveat)
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
    MMS = sim.MMS

    nn_index = 1

    Te_index = nvariables - 2
    ne_index = nvariables - 1
    ϕ_index  = nvariables

    mi = sim.propellant.m

    for (i, z) in enumerate(z_cell)
        if MMS
            U[nn_index, i] = sim.initial_nn_mms(z)
        else
            U[nn_index, i] = 1.5e19 * mi
        end

        ni = sim.initial_ni(z)
        Te = sim.initial_Te(z)
        ϕ = sim.initial_ϕ(z)
        Ti = sim.initial_Ti(z)

        E = sim.propellant.cv * Ti
      

        # ions initialized with equal densities, same velocity as neutrals
        for j in fluid_ranges[2:end]
            n_index = j[1]
            nu_index = j[2]
            nE_index = j[3]
            U[n_index, i] = ni * mi
            U[nu_index, i] = ni * un * mi
            U[nE_index, i] = ni * E * mi
        end
        U[Te_index, i] = Te
        @views U[ne_index, i] = electron_density(U[:, i], fluid_ranges)
        U[ϕ_index, i]  = ϕ
    end
    
    return U
end

function initial_condition_v2!(U, z_cell, ic, fluid_ranges) #start with only one fluid, can be extended to multiple fluids later
    #can extend later to more 
    #also not using inlet_neutral_density for now
    nvariables = size(U, 1)
    #nn = inlet_neutral_density(sim)
    mi = HallThruster.Air.m
    
    for (i, z) in enumerate(z_cell)
        # ions initialized with equal densities, same velocity as neutrals
        for (m, j) in enumerate(fluid_ranges)
            n = ic[1](z)
            u = ic[2](z)
            T = ic[3](z)
            E = HallThruster.Air.cv * T
            n_index = j[1]
            U[n_index, i] = n * mi
            if j[2] !== nothing
                nu_index = j[2]
                U[nu_index, i] = n * u * mi
            end
            if j[3] !== nothing
                nE_index = j[3]
                U[nE_index, i] = n * E * mi
            end
        end
    end
    return U
end