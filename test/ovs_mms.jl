using HallThruster, StaticArrays, Symbolics, DifferentialEquations, Statistics, Plots

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

function perform_OVS(; MMS_CONSTS, fluxfn, reconstruct)

    source = mms!
    
    #need to somehow get the functions from runtests.jl to here
    #think about how to integrate ICs and BCs such that no need to define quantities twice, often BCs will be similar to ICs.
    #need information about fluid in BCs for certain quantities as well. 

    function IC!(U, z, fluids, L)
        gas1 = fluids[1].species.element
        gas2 = fluids[2].species.element
        #gas3 = fluids[3].species.element

        ρ1 = 2000.0
        ρ2 = 3000.0
        u2 = 300.0
        ρ3 = 5.0
        u3 = 200.0
        T3 = 300.0

        U .= [ρ1, ρ2, ρ2*u2] #, ρ3, ρ3*u3, ρ3*T3*gas3.cv]
        return U
    end

    ρ1 = 2000.0
    ρ2 = 2000.0
    u2 = 300.0
    ρ3 = 5.0
    u3 = 200.0
    T3 = 300.0

    left_state = [ρ1, ρ2, ρ2*u2] #, ρ3, ρ3*u3, ρ3*T3*HallThruster.Xenon.cv]
    right_state = [ρ1, ρ2+1000, (ρ2+1000)*(u2+100)] #, ρ3, ρ3*u3, ρ3*T3*HallThruster.Xenon.cv]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Dirichlet(right_state))

    simulation = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(HallThruster.SPT_100, MMS_CONSTS.n_cells_start), 
    boundary_conditions = BCs,
    scheme = HallThruster.HyperbolicScheme(fluxfn, HallThruster.minmod, reconstruct),
    initial_condition = IC!, 
    source_term! = source, 
    fluids = [HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.ContinuityOnly(MMS_CONSTS.u_constant, MMS_CONSTS.T_constant));
    HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.IsothermalEuler(MMS_CONSTS.T_constant))], #;
    #HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())], 
    end_time = MMS_CONSTS.max_end_time, 
    saveat = [MMS_CONSTS.max_end_time],
    timestepcontrol = (1e-6, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
    callback = DiffEqCallbacks.TerminateSteadyState(1e-8, 1e-6, DiffEqCallbacks.allDerivPass) 
    )

    #generate different number of cells while keeping CFL number constant over refinements chosen
    refinements = MMS_CONSTS.refinements
    results = Array{Result, 1}(undef, refinements)
    simulations_mms = Array{HallThruster.MultiFluidSimulation, 1}(undef, refinements)
    n_cells = MMS_CONSTS.n_cells_start
    _, __, fluid_ranges, ___ = HallThruster.configure_simulation(simulation)

    for refinement in 1:refinements #set timestep and iterate number of cells
        simulation.grid = HallThruster.generate_grid(HallThruster.SPT_100, n_cells)
        simulation.timestepcontrol = (MMS_CONSTS.CFL*MMS_CONSTS.L/(MMS_CONSTS.u_constant*n_cells), false)
        simulations_mms[refinement] = deepcopy(simulation)
        n_cells = n_cells*2
    end
    
    
    for refinement in 1:refinements
        sol = HallThruster.run_simulation(simulations_mms[refinement])
        z_cells = simulations_mms[refinement].grid.cell_centers
        u_exa = Array{Union{Nothing, Float64}}(nothing, length(sol.u[1][:, 1]), length(z_cells))
        for (i, z_cell) in enumerate(z_cells)
            u_exa[1, i] = nn_manufactured_f(z_cell, MMS_CONSTS)
            for (j, index) in enumerate(fluid_ranges[2:end])
                u_exa[index[1], i] = ni_manufactured_f(z_cell, MMS_CONSTS)
                u_exa[index[2], i] = ni_manufactured_f(z_cell, MMS_CONSTS) * ui_manufactured_f(z_cell, MMS_CONSTS)
            end
        end
        
        error = abs.(u_exa - sol.u[1])
        results[refinement] = Result(sol, simulations_mms[refinement].timestepcontrol[1], z_cells, simulations_mms[refinement].grid.ncells, u_exa, error, [maximum(error[i, :]) for i in 1:size(error)[1]], [Statistics.mean(error[i, :]) for i in 1:size(error)[1]])
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
    L_1 = Array{Union{Nothing, Float64}}(nothing, 3)
    L_inf = Array{Union{Nothing, Float64}}(nothing, 3)
    for j in 1:3  ###need to rewrite this without ncharge terms
        L_1[j] = compute_slope([results[i].ncells for i in 1:length(results)], [results[i].L_1[j] for i in 1:length(results)])
        L_inf[j] = compute_slope([results[i].ncells for i in 1:length(results)], [results[i].L_inf[j] for i in 1:length(results)])
    end
    return L_1, L_inf
end

