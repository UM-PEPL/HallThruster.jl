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

function perform_OVS(; MMS_CONSTS, fluxfn, reconstruct)

    #create a template definition of source! function somewhere
    function source!(Q, U, params, ϕ, Tev, i)
        mms!(Q, [params.z_cell[i]])
    end
    
    function IC!(U, z, fluids, L)
        gas1 = fluids[1].species.element

        ρ1 = 3000.0
        u1 = MMS_CONSTS.u_constant
        T1 = MMS_CONSTS.T_constant
        E = MMS_CONSTS.fluid.cv*T1 + 0.5*u1*u1

        U .= [ρ1, ρ1, ρ1*u1] #[ρ1, ρ1*u1, ρ1*E]f
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

    simulation = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(HallThruster.SPT_100, MMS_CONSTS.n_cells_start), 
    boundary_conditions = BCs,
    scheme = HallThruster.HyperbolicScheme(fluxfn, HallThruster.minmod, reconstruct),
    initial_condition = IC!, 
    source_term! = source!, 
    fluids = [HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.ContinuityOnly(MMS_CONSTS.u_constant, MMS_CONSTS.T_constant));
    HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.IsothermalEuler(MMS_CONSTS.T_constant))],
    #[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())], 
    end_time = MMS_CONSTS.max_end_time, 
    saveat = [MMS_CONSTS.max_end_time],
    timestepcontrol = (1e-6, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
    callback = DiffEqCallbacks.TerminateSteadyState(1e-8, 1e-6, DiffEqCallbacks.allDerivPass) 
    )

    #generate different number of cells while keeping CFL number constant over refinements chosen
    refinements = MMS_CONSTS.refinements
    results = Array{Result, 1}(undef, refinements)
    n_cells = MMS_CONSTS.n_cells_start

    for refinement in 1:refinements
        simulation.grid = HallThruster.generate_grid(HallThruster.SPT_100, n_cells)
        simulation.timestepcontrol = (MMS_CONSTS.CFL*MMS_CONSTS.L/(MMS_CONSTS.u_constant*n_cells), false)
        sol = HallThruster.run_simulation(simulation)
        z_cells = simulation.grid.cell_centers
        u_exa = Array{Union{Nothing, Float64}}(nothing, length(sol.u[1][:, 1]), length(z_cells))
        for (i, z_cell) in enumerate(z_cells)
            u_exa[:, i] = mms_conservative([z_cell, 0])
        end
        error = abs.(u_exa - sol.u[1])
        results[refinement] = Result(sol, simulation.timestepcontrol[1], z_cells, simulation.grid.ncells, u_exa, error, [maximum(error[i, :]) for i in 1:size(error)[1]], [Statistics.mean(error[i, :]) for i in 1:size(error)[1]])
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
    L_1 = Array{Union{Nothing, Float64}}(nothing, nfluids)
    L_inf = Array{Union{Nothing, Float64}}(nothing, nfluids)
    for j in 1:nfluids
        L_1[j] = compute_slope([results[i].ncells for i in 1:length(results)], [results[i].L_1[j] for i in 1:length(results)])
        L_inf[j] = compute_slope([results[i].ncells for i in 1:length(results)], [results[i].L_inf[j] for i in 1:length(results)])
    end
    return L_1, L_inf
end

