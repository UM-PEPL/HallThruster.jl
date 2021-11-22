
using Test, HallThruster, Plots

function source!(Q, U, params, ϕ, i)
    HallThruster.apply_reactions!(Q, U, params, i)
    HallThruster.apply_ion_acceleration!(Q, U, params, ϕ, i)
    return Q
end

function IC!(U, z, fluids, L)
    ρ1 = 1.0
    ρ2 = 0.01
    u1 = 300.0
    U[1] = ρ1
    U[2] = ρ2
    U .= SA[ρ1, ρ2, ρ2*u1] #[ρ1, ρ1*u1, ρ1*E]
    return U
end


function run(end_time = 0.0002)
    fluid = HallThruster.Xenon
    timestep = 1e-8

    ρ1 = 1.0
    ρ2 = 0.01
    u1 = 300.0
    T1 = 300.0

    left_state = [ρ1, ρ2, ρ2*u1] # [ρ1, ρ1*u1, ρ1*E]
    right_state = [ρ1, ρ2, ρ2*(u1+0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Neumann())

    sim = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(HallThruster.SPT_100, 100),
        boundary_conditions = BCs,
        scheme = HallThruster.HyperbolicScheme(HallThruster.HLLE!, identity, false),
        initial_condition = IC!, source_term! = source!,
        fluids = [HallThruster.Fluid(HallThruster.Species(fluid, 0), HallThruster.ContinuityOnly(300.0, 300.0));
        HallThruster.Fluid(HallThruster.Species(fluid, 1), HallThruster.IsothermalEuler(300.0))],
        #[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())],
        end_time = end_time, #0.0002
        saveat = [end_time],
        timestepcontrol = (timestep, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
        callback = nothing
        )

    @time sol = HallThruster.run_simulation(sim)

    #extract potential at the end, just for now, make proper later ##############################################################
    species, fluids, fluid_ranges, species_range_dict = HallThruster.configure_simulation(sim)
    grid = sim.grid

    U, cache = HallThruster.allocate_arrays(sim)

    HallThruster.initial_condition!(U, grid.cell_centers, sim.initial_condition, fluid_ranges, fluids)

    scheme = sim.scheme
    source_term! = sim.source_term!
    timestep = sim.timestepcontrol[1]
    adaptive = sim.timestepcontrol[2]
    tspan = (0., sim.end_time)

    reactions = HallThruster.load_ionization_reactions(species)
    BCs = sim.boundary_conditions

    params = (;
        cache,
        fluids,
        fluid_ranges,
        species_range_dict,
        z_cell = grid.cell_centers,
        z_edge = grid.edges,
        cell_volume = grid.cell_volume,
        source_term!,
        reactions,
        scheme,
        BCs,
        dt = timestep
    )

    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    reactions, species_range_dict = params.reactions, params.species_range_dict

    F, UL, UR, Q, A, b, ϕ, Tev = params.cache

    z_cell, z_edge, cell_volume = params.z_cell, params.z_edge, params.cell_volume
    scheme = params.scheme
    source_term! = params.source_term!

    nvariables = size(U, 1)
    ncells = size(U, 2) - 2

    #make U last timestep
    #U = sol.u[1]
    Tev .= 5.0

    A .= 0.0
    b .= 0.0
    HallThruster.set_up_potential_equation!(U, A, b, Tev, params)
    ϕ = A\b

end