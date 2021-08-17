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

function run_simulation(sim)

    fluids, fluid_ranges, species_range_dict = configure_simulation(sim)
    z_cell, z_edge = generate_grid(sim.geometry, sim.ncells)

    U, cache = allocate_arrays(sim)

    IC = initial_condition!(U, z_cell, sim, fluid_ranges)

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
