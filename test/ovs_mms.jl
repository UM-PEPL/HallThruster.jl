using HallThruster, StaticArrays, Symbolics, OrdinaryDiffEq, DiffEqBase, Statistics, Plots, DiffEqCallbacks

mutable struct Result
    solution #::Vector{Matrix{Float64}}
    timestep
    z_cells #::Vector{Float64}
    ncells
    u_exa
    errors
    L_inf
    L_1
end

"""
    perform_OVS(; MMS_CONSTS::NamedTuple, fluxfn::Function, reconstruct::Bool)

function to perform OVS for the fluid equations using definition of manufactured
solution in runtests.jl
"""

function perform_OVS(; MMS_CONSTS, fluxfn, reconstruct)

    #create a template definition of source! function somewhere
    function source!(Q, U, params, i)
        mms!(@views(Q[1:3]), [params.z_cell[i]])
        Q[4] = 0
    end

    function source_potential!(b, U, s_consts)
        HallThruster.potential_source_term!(b, U, s_consts)
        #HallThruster.OVS_potential_source_term!(b, s_consts)
    end
    
    function boundary_potential!(A, b, U, bc_consts)
        ϕ_L = 400.0
        ϕ_R = 0.0
        HallThruster.boundary_conditions_potential!(A, b, U, bc_consts, ϕ_L, ϕ_R)
        #HallThruster.OVS_boundary_conditions_potential!((A, b, U, bc_consts, ϕ_L, ϕ_R)
    end
    
    function IC!(U, z, fluids, L)
        gas1 = fluids[1].species.element

        ρ1 = 3000.0
        u1 = 300.0
        #u1 = MMS_CONSTS.u_constant

        U .= [ρ1, ρ1, ρ1*u1, 0.0] #[ρ1, ρ1*u1, ρ1*E]f
        return U
    end

    ρ1 = MMS_CONSTS.n0 + MMS_CONSTS.nx
    u1 = MMS_CONSTS.u_constant
    T1 = MMS_CONSTS.T_constant
    E = MMS_CONSTS.fluid.cv*T1 + 0.5*u1*u1
    ER = MMS_CONSTS.fluid.cv*(T1+100.0) + 0.5*(u1+0.0)*(u1+0.0)

    left_state = [ρ1, ρ1, ρ1*u1] # [ρ1, ρ1*u1, ρ1*E] 
    right_state = [ρ1, ρ1, ρ1*(u1+0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Neumann())

    n_save = 1
    saveat = if n_save == 1
        [MMS_CONSTS.max_end_time]
    else
        LinRange(0.0, MMS_CONSTS.max_end_time, n_save) |> collect
    end

    simulation = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(HallThruster.SPT_100, MMS_CONSTS.n_cells_start), 
    boundary_conditions = BCs,
    scheme = HallThruster.HyperbolicScheme(fluxfn, HallThruster.minmod, reconstruct),
    initial_condition = IC!, 
    source_term! = source!,
    source_potential! = source_potential!,
    boundary_potential! = boundary_potential!, 
    fluids = [HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.ContinuityOnly(MMS_CONSTS.u_constant, MMS_CONSTS.T_constant));
    HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.IsothermalEuler(MMS_CONSTS.T_constant))],
    #[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())], 
    end_time = MMS_CONSTS.max_end_time, 
    saveat = saveat,
    timestepcontrol = (1e-6, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
    callback = DiffEqCallbacks.TerminateSteadyState(1e-8, 1e-6, DiffEqCallbacks.allDerivPass), 
    solve_energy = false
    )

    #generate different number of cells while keeping CFL number constant over refinements chosen
    refinements = MMS_CONSTS.refinements
    results = Array{Result, 1}(undef, refinements)
    n_cells = MMS_CONSTS.n_cells_start
    _, fluids, fluid_ranges, __ = HallThruster.configure_simulation(simulation)
    lf = fluid_ranges[end][end]

    for refinement in 1:refinements
        simulation.grid = HallThruster.generate_grid(HallThruster.SPT_100, n_cells)
        simulation.timestepcontrol = (MMS_CONSTS.CFL*MMS_CONSTS.L/(MMS_CONSTS.u_constant*n_cells), false)
        sol = HallThruster.run_simulation(simulation)
        z_cells = simulation.grid.cell_centers
        u_exa = Array{Union{Nothing, Float64}}(nothing, length(sol.u[1][1:lf, 1]), length(z_cells))
        for (i, z_cell) in enumerate(z_cells)
            u_exa[:, i] = mms_conservative([z_cell, 0])
        end
        error = abs.(u_exa - sol.u[end][1:lf, :])./u_exa[:, 1]
        results[refinement] = Result(sol, simulation.timestepcontrol[1], z_cells, simulation.grid.ncells, u_exa, error, [maximum(error[i, :]) for i in 1:size(error)[1]], [Statistics.mean(error[i, :]) for i in 1:size(error)[1]])
        n_cells = n_cells*2
    end
    return results
end

function perform_OVS_potential(; MMS_CONSTS, fluxfn, reconstruct)

    #need to loop through grid refinements
    refinements = MMS_CONSTS.refinements
    results = Array{Result, 1}(undef, refinements)
    n_cells = MMS_CONSTS.n_cells_start

    for refinement in 1:refinements
        U, params = set_up_params_U(; MMS_CONSTS, fluxfn, reconstruct, ncells = n_cells)
        HallThruster.solve_potential!(U, params)
        z_cells = params.z_cell
        ϕ = U[params.index.ϕ, :]
        u_exa = Array{Union{Nothing, Float64}}(nothing, length(ϕ)-2)
        for i in 1:length(ϕ)-2
            u_exa[i] = 25000.0*z_cells[i+1]^2 - 3250.0*z_cells[i+1] + 100.0
        end
        error = abs.(u_exa - ϕ[2:end-1])
        #timestep irrelevant, adapt ncells to amount of cells for 
        results[refinement] = Result(ϕ, 1.0, z_cells, n_cells-2, u_exa, error, [maximum(error[i, :]) for i in 1:size(error)[1]], [Statistics.mean(error[i, :]) for i in 1:size(error)[1]])
        n_cells = n_cells*2
    end
    return results
end

function compute_slope(ncells, errors)
    p = Array{Union{Nothing, Float64}}(nothing, length(ncells)-2)
    for i in 1:length(ncells)-2
        p[i] = log(abs(errors[i+2]-errors[i+1])/abs(errors[i+1]-errors[i]))/log(0.5)
    end 
    return Statistics.mean(p)
end

function evaluate_slope(results, MMS_CONSTS)
    nfluids = size(results[1].u_exa)[1]
    if nfluids > MMS_CONSTS.refinements #- 1 #for potential
        nfluids = 1
    end
    L_1 = Array{Union{Nothing, Float64}}(nothing, nfluids)
    L_inf = Array{Union{Nothing, Float64}}(nothing, nfluids)
    for j in 1:nfluids
        L_1[j] = compute_slope([results[i].ncells for i in 1:length(results)], [results[i].L_1[j] for i in 1:length(results)])
        L_inf[j] = compute_slope([results[i].ncells for i in 1:length(results)], [results[i].L_inf[j] for i in 1:length(results)])
    end
    return L_1, L_inf
end

"""
    set_up_params_U(; MMS_CONSTS, fluxfn, reconstruct)

just a helper function to generate U and params for use in potential OVS. Not useful for anything else.
"""

function set_up_params_U(; MMS_CONSTS, fluxfn, reconstruct, ncells) #does not really matter, just to generate params and U
    
    function source!(Q, U, params, i)
    end

    function source_potential!(b, U, s_consts)
        #HallThruster.potential_source_term!(b, U, s_consts)
        HallThruster.OVS_potential_source_term!(b, s_consts)
    end
    
    function boundary_potential!(A, b, U, bc_consts)
        ϕ_L = 100.0
        ϕ_R = 0.0
        #HallThruster.boundary_conditions_potential!(A, b, U, bc_consts, ϕ_L, ϕ_R)
        HallThruster.OVS_boundary_conditions_potential!(A, b, U, bc_consts, ϕ_L, ϕ_R)
    end

    function IC!(U, z, fluids, L)
        gas1 = fluids[1].species.element

        ρ1 = 3000.0
        u1 = MMS_CONSTS.u_constant
        T1 = MMS_CONSTS.T_constant
        E = MMS_CONSTS.fluid.cv*T1 + 0.5*u1*u1

        U .= [ρ1, ρ1, ρ1*u1, 0.0] #[ρ1, ρ1*u1, ρ1*E]
        return U
    end

    ρ1 = 3000.0
    u1 = MMS_CONSTS.u_constant
    T1 = MMS_CONSTS.T_constant
    E = MMS_CONSTS.fluid.cv*T1 + 0.5*u1*u1
    ER = MMS_CONSTS.fluid.cv*(T1+100.0) + 0.5*(u1+0.0)*(u1+0.0)

    left_state = [ρ1, ρ1, ρ1*u1] # [ρ1, ρ1*u1, ρ1*E] 
    right_state = [ρ1, ρ1, ρ1*(u1+0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Dirichlet(right_state))

    simulation = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(HallThruster.SPT_100, ncells), 
    boundary_conditions = BCs,
    scheme = HallThruster.HyperbolicScheme(fluxfn, HallThruster.minmod, reconstruct),
    initial_condition = IC!, 
    source_term! = source!,
    source_potential! = source_potential!,
    boundary_potential! = boundary_potential!, 
    fluids = [HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.ContinuityOnly(MMS_CONSTS.u_constant, MMS_CONSTS.T_constant));
    HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.IsothermalEuler(MMS_CONSTS.T_constant))],
    #[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())], 
    end_time = MMS_CONSTS.max_end_time, 
    saveat = [MMS_CONSTS.max_end_time],
    timestepcontrol = (1e-6, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
    callback = nothing,
    solve_energy = false
    )

    sim = simulation

    species, fluids, fluid_ranges, species_range_dict = HallThruster.configure_simulation(sim)
    grid = sim.grid

    U, cache = HallThruster.allocate_arrays(sim)

    lf = fluid_ranges[end][end]
    index = (;lf = lf, nϵ = lf+1, Tev = lf+2, ne = lf+3, pe = lf+4, ϕ = lf+5)

    HallThruster.initial_condition!(@views(U[1:index.nϵ, :]), grid.cell_centers, sim.initial_condition, fluids)

    scheme = sim.scheme
    source_term! = sim.source_term!
    timestep = sim.timestepcontrol[1]
    adaptive = sim.timestepcontrol[2]
    tspan = (0.0, sim.end_time)

    reactions = HallThruster.load_ionization_reactions(species)
    landmark = HallThruster.load_landmark()

    BCs = sim.boundary_conditions

    HallThruster.precompute_bfield!(cache.B, grid.cell_centers)

    params = (; index, cache, fluids, fluid_ranges, species_range_dict, z_cell=grid.cell_centers,
              z_edge=grid.edges, cell_volume=grid.cell_volume, source_term!, reactions,
              scheme, BCs, dt=timestep, source_potential! = sim.source_potential!, 
              boundary_potential! = sim.boundary_potential!, landmark, solve_energy = sim.solve_energy)

    return U, params
end

function perform_OVS_elecenergy(; MMS_CONSTS, fluxfn, reconstruct)

    #create a template definition of source! function somewhere
    function source!(Q, U, params, i)
        mms!(@views(Q[1:4]), [params.z_cell[i]])
        #Q[4] = 0.0
    end

    function source_potential!(b, U, s_consts)
        HallThruster.potential_source_term!(b, U, s_consts)
        #HallThruster.OVS_potential_source_term!(b, s_consts)
    end
    
    function boundary_potential!(A, b, U, bc_consts)
        ϕ_L = 400.0
        ϕ_R = 0.0
        HallThruster.boundary_conditions_potential!(A, b, U, bc_consts, ϕ_L, ϕ_R)
        #HallThruster.OVS_boundary_conditions_potential!((A, b, U, bc_consts, ϕ_L, ϕ_R)
    end
    
    function IC!(U, z, fluids, L)
        ρ2 = MMS_CONSTS.n0 + MMS_CONSTS.nx
        u1 = MMS_CONSTS.u0
        ρ1 = MMS_CONSTS.n0 + MMS_CONSTS.nx
        #Tev = 3.0
        #Tev = MMS_CONSTS.Tev_elec_max*exp(-(2 * (z - MMS_CONSTS.L/2) / 0.033)^2)
        Tev = MMS_CONSTS.Tev0 + MMS_CONSTS.Tev_elec_max*sin(π * z / (MMS_CONSTS.L))
        ne = (MMS_CONSTS.n0 + MMS_CONSTS.nx) / fluids[1].species.element.m
        U .= SA[ρ1, ρ2, ρ2*u1, 3/2*ne*Tev*HallThruster.kB]
        return U
    end

    u1 = MMS_CONSTS.u0
    ρ1 = MMS_CONSTS.n0 + MMS_CONSTS.nx
    T1 = MMS_CONSTS.T_constant
    E = MMS_CONSTS.fluid.cv*T1 + 0.5*u1*u1
    ER = MMS_CONSTS.fluid.cv*(T1+100.0) + 0.5*(u1+0.0)*(u1+0.0)

    left_state = [ρ1, ρ1, ρ1*u1] # [ρ1, ρ1*u1, ρ1*E] 
    right_state = [ρ1, ρ1, ρ1*(u1+0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Dirichlet(right_state))

    n_save = 100
    saveat = if n_save == 1
        [MMS_CONSTS.max_end_time]
    else
        LinRange(0.0, MMS_CONSTS.max_end_time, n_save) |> collect
    end

    condition(u,t,integrator) = t < 1
    function affect!(integrator)
        U, params = integrator.u, integrator.p
        
        fluids, fluid_ranges = params.fluids, params.fluid_ranges
        index = params.index

        B = params.cache.B

        z_cell, z_edge, cell_volume = params.z_cell, params.z_edge, params.cell_volume
        ncells = size(U, 2) - 2

        ####################################################################
        #PREPROCESS
        #calculate useful quantities relevant for potential, electron energy and fluid solve
        L_ch = 0.025
        fluid = fluids[1].species.element

        @inbounds for i in 1:(ncells + 2)
            #update electron temperature from energy using old density
            if params.solve_energy
                #U[index.Tev, i] = max(1, U[index.nϵ, i]/3*2/U[index.ne, i])
            end
            U[index.ne, i] = max(1e-10, HallThruster.electron_density(@view(U[:, i]), fluid_ranges) / fluid.m)
            U[index.pe, i] = HallThruster.electron_pressure(U[index.ne, i], U[index.Tev, i]) #this would be real electron pressure, ie next step use for previous in energy convection update
            #U[index.pe, i] = U[index.nϵ, i]/3*2*HallThruster.e #if using the same for pe and ne, might solve some instabilities
            U[index.grad_ϕ, i] = HallThruster.first_deriv_central_diff(U[index.ϕ, :], params.z_cell, i)
            U[index.ue, i] = -0.0001 #HallThruster.electron_velocity(U, params, i) #for first try, set equal to 2000
            params.cache.νan[i] = HallThruster.get_v_an(z_cell[i], B[i], L_ch)
            params.cache.νc[i] = HallThruster.get_v_c(U[index.Tev, i], U[1, i]/fluid.m , U[index.ne, i], fluid.m)
            params.cache.μ[i] = HallThruster.cf_electron_transport(params.cache.νan[i], params.cache.νc[i], B[i])
        end
        
        #POTENTIAL #########################################################
        HallThruster.solve_potential!(U, params)
    end

    cb_update = DiscreteCallback(condition, affect!, save_positions=(false,false))
    cb_convergence = DiffEqCallbacks.TerminateSteadyState(1e-8, 1e-6, DiffEqCallbacks.allDerivPass)
    cb = CallbackSet(cb_update, cb_convergence)

    simulation = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(HallThruster.SPT_100, MMS_CONSTS.n_cells_start), 
    boundary_conditions = BCs,
    scheme = HallThruster.HyperbolicScheme(fluxfn, HallThruster.minmod, reconstruct),
    initial_condition = IC!, 
    source_term! = source!,
    source_potential! = source_potential!,
    boundary_potential! = boundary_potential!, 
    fluids = [HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.ContinuityOnly(MMS_CONSTS.u_constant, MMS_CONSTS.T_constant));
    HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.IsothermalEuler(MMS_CONSTS.T_constant))],
    #[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())], 
    end_time = MMS_CONSTS.max_end_time, 
    saveat = saveat,
    timestepcontrol = (1e-6, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
    callback = cb, 
    solve_energy = true
    )

    #generate different number of cells while keeping CFL number constant over refinements chosen
    refinements = MMS_CONSTS.refinements
    results = Array{Result, 1}(undef, refinements)
    n_cells = MMS_CONSTS.n_cells_start
    _, fluids, fluid_ranges, __ = HallThruster.configure_simulation(simulation)
    lf = fluid_ranges[end][end]

    for refinement in 1:refinements
        simulation.grid = HallThruster.generate_grid(HallThruster.SPT_100, n_cells)
        simulation.timestepcontrol = (MMS_CONSTS.CFL*MMS_CONSTS.L/(MMS_CONSTS.u_constant*n_cells), false)
        sol = HallThruster.run_simulation(simulation)
        z_cells = simulation.grid.cell_centers
        u_exa = Array{Union{Nothing, Float64}}(nothing, length(sol.u[1][1:lf+1, 1]), length(z_cells)) #lf+1 for the electron energy
        for (i, z_cell) in enumerate(z_cells)
            u_exa[:, i] = mms_conservative([z_cell, 0])
        end
        error = abs.(u_exa - sol.u[end][1:lf+1, :])./u_exa[:, 1]
        results[refinement] = Result(sol, simulation.timestepcontrol[1], z_cells, simulation.grid.ncells, u_exa, error, [maximum(error[i, :]) for i in 1:size(error)[1]], [Statistics.mean(error[i, :]) for i in 1:size(error)[1]])
        n_cells = n_cells*2
    end
    return results
end

#eventually do ovs of whole system together without having to do it separately for each equation. 
