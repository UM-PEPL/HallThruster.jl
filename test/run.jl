
using Test, HallThruster, Plots, StaticArrays

function source!(Q, U, params, ϕ, Tev, i)
    HallThruster.apply_reactions!(Q, U, params, Tev, i)
    HallThruster.apply_ion_acceleration!(Q, U, params, ϕ, i)
    return Q
end

function IC!(U, z, fluids, L)
    ρ1 = 2.1801715574645586e-6
    ρ2 = ρ1 * exp(-((z - L)/0.033)^2)
    u1 = 300.0
    U[1] = ρ1
    U[2] = ρ2
    U .= SA[ρ1, ρ2, ρ2*u1] #[ρ1, ρ1*u1, ρ1*E]
    return U
end

function run(end_time = 0.0002, n_save = 2)
    fluid = HallThruster.Xenon
    timestep = 0.9e-8

    ρ1 = 2.1801715574645586e-6
    ρ2 = 2.1801715574645586e-6
    u1 = 300.0
    T1 = 1000.0

    left_state = [ρ1, ρ2, ρ2*u1] # [ρ1, ρ1*u1, ρ1*E]
    right_state = [ρ1, ρ2, ρ2*(u1+0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Neumann())

    sim = HallThruster.MultiFluidSimulation(
        grid = HallThruster.generate_grid(HallThruster.SPT_100, 100),
        boundary_conditions = BCs,
        scheme = HallThruster.HyperbolicScheme(HallThruster.HLLE!, identity, false),
        initial_condition = IC!, source_term! = source!,
        fluids = [HallThruster.Fluid(HallThruster.Species(fluid, 0), HallThruster.ContinuityOnly(300.0, 300.0));
        HallThruster.Fluid(HallThruster.Species(fluid, 1), HallThruster.IsothermalEuler(300.0))],
        #[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())],
        end_time = end_time, #0.0002
        saveat = if n_save == 1
                [end_time]
            else
                LinRange(0.0, end_time, n_save) |> collect
            end,
        timestepcontrol = (timestep, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
        callback = nothing
    )

    @time sol = HallThruster.run_simulation(sim)

    p = plot()#plot(sol.u[end][1, :], yaxis = :log)
    plot!(p, sol.u[end][3, :] ./ sol.u[end][2, :])

    display(p)

    #@show fieldnames(typeof(sol))
    return sol

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
    Tev .= 2.0

    A .= 0.0
    b .= 0.0
    HallThruster.set_up_potential_equation!(U, A, b, Tev, params)
    ϕ = A\b

    plot(z_cell, ϕ)

end

function animate_solution(sol)
    mi = HallThruster.Xenon.m
    @gif for (u, t) in zip(sol.u, sol.t)
        p = plot(ylims = (1e13, 1e20))
        plot!(p, u[1, :]/mi, yaxis = :log)
        plot!(p, u[2, :]/mi)
    end
end