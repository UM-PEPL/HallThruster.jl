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
    
    simulation = (
    ncells = 100,
    propellant = HallThruster.Xenon,
    ncharge = MMS_CONSTS.ncharge,
    MMS = true,
    mms! = mms!, 
    cb = DiffEqCallbacks.TerminateSteadyState(1e-8, 1e-6, DiffEqCallbacks.allDerivPass),  #abstol, reltol, test
    geometry = HallThruster.SPT_100,
    neutral_temperature = 500.,
    neutral_velocity = MMS_CONSTS.un,
    ion_temperature = MMS_CONSTS.ion_temperature,
    initial_nn_mms = nn_mms_func, 
    initial_Te = Te_func,
    initial_ϕ = ϕ_func,
    initial_ni = ni_func,
    solve_Te = false,
    solve_ne = false,
    inlet_mdot = 5e-6,
    saveat = [MMS_CONSTS.max_end_time],
    tspan = (0., MMS_CONSTS.max_end_time),
    dt = 200e-8, #5e-8
    scheme = (
        flux_function = fluxfn,
        limiter = identity,
        reconstruct = reconstruct
        ),
    )

    refinements = MMS_CONSTS.refinements
    results = Array{Result, 1}(undef, refinements)
    simulations_mms = Array{NamedTuple, 1}(undef, refinements)
    n_cells = MMS_CONSTS.n_cells_start
    _, __, fluid_ranges, ___ = HallThruster.configure_simulation(simulation)

    for refinement in 1:refinements
        simulations_mms[refinement] = merge(simulation, (; ncells = n_cells, dt = MMS_CONSTS.CFL*MMS_CONSTS.L/(MMS_CONSTS.un*n_cells)))
        n_cells = n_cells*2
    end

    for refinement in 1:refinements
        sol = HallThruster.run_simulation(simulations_mms[refinement])
        z_cells = HallThruster.generate_grid(simulation.geometry, simulations_mms[refinement].ncells)[1]
        u_exa = Array{Union{Nothing, Float64}}(nothing, length(sol.u[1][:, 1]), length(z_cells))
        for (i, z_cell) in enumerate(z_cells)
            u_exa[1, i] = nn_manufactured_f(z_cell, MMS_CONSTS)
            for (j, index) in enumerate(fluid_ranges[2:end])
                u_exa[index[1], i] = ni_manufactured_f(z_cell, MMS_CONSTS)
                u_exa[index[2], i] = ni_manufactured_f(z_cell, MMS_CONSTS) * ui_manufactured_f(z_cell, MMS_CONSTS)
            end
        end
        u_exa[end-2:end, :] .= 0.0
        error = abs.(u_exa - sol.u[1])
        results[refinement] = Result(sol, simulations_mms[refinement].dt, z_cells, simulations_mms[refinement].ncells, u_exa, error, [maximum(error[i, :]) for i in 1:size(error)[1]], [Statistics.mean(error[i, :]) for i in 1:size(error)[1]])
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
    L_1 = Array{Union{Nothing, Float64}}(nothing, MMS_CONSTS.ncharge*2+1)
    L_inf = Array{Union{Nothing, Float64}}(nothing, MMS_CONSTS.ncharge*2+1)
    for j in 1:MMS_CONSTS.ncharge*2+1
        L_1[j] = compute_slope([results[i].ncells for i in 1:length(results)], [results[i].L_1[j] for i in 1:length(results)])
        L_inf[j] = compute_slope([results[i].ncells for i in 1:length(results)], [results[i].L_inf[j] for i in 1:length(results)])
    end
    return L_1, L_inf
end

