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
    initial_condition::IC
    boundary_conditions::Tuple{B1, B2}   # Tuple of left and right boundary conditions, subject to the approval of PR #10 
    end_time::Float64    # How long to simulate
    scheme:: HyperbolicScheme{F, L} # Flux, Limiter
    source_term!::S  # Source term function. This can include reactons, electric field, and MMS terms
    saveat::Vector{Float64} #when to save
end

get_species(sim) = [Species(sim.propellant, i) for i in 0:sim.ncharge]

function configure_simulation(sim)
    fluids = sim.fluids
    species = [fluids[i].species for i in 1:length(fluids)]
    fluid_ranges = ranges(fluids)
    species_range_dict = Dict(
        fluid.species => fluid_range for (fluid, fluid_range) in zip(fluids, fluid_ranges)
    )

    return species, fluids, fluid_ranges, species_range_dict
end

function allocate_arrays(sim) #rewrite allocate arrays as function of set of equations, either 1, 2 or 3
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

function update!(dU, U, params, t)
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
    grid = sim.grid

    U, cache = allocate_arrays(sim)

    initial_condition!(U, grid.cell_centers, sim.initial_condition, fluid_ranges, fluids)

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

    prob = ODEProblem{true}(update!, U, tspan, params)
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

function initial_condition!(U, z_cell, IC!, fluid_ranges, fluids)
    #can extend later to more 
    #also not using inlet_neutral_density for now
    #nn = inlet_neutral_density(sim)
    
    for (i, z) in enumerate(z_cell)
        @views IC!(U[:, i], z, fluids, z_cell[end])
    end
    return U
end