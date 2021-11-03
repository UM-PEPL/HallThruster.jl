using Test, Documenter, HallThruster, StaticArrays, BenchmarkTools, Symbolics, DifferentialEquations, Statistics, Plots

Te_func = z -> 30 * exp(-(2(z - HallThruster.SPT_100.channel_length) / 0.033)^2)
ϕ_func = z -> 300 * (1 - 1/(1 + exp(-1000 * (z - HallThruster.SPT_100.channel_length))))
ni_func = z -> 2000 #1e6
nn_mms_func = z -> 2000

end_time = 20e-5 #30e-5

const MMS_CONSTS = (
    CFL = 0.99, 
    n_cells_start = 10,
    refinements = 7,
    n_waves = 2.0,
    un = 300.0, 
    L = HallThruster.SPT_100.domain[2]-HallThruster.SPT_100.domain[1],
    ion_temperature = 0.0,
    nn0 = 1000.0,
    nnx = 1000.0,
    ni0 = 2000.0,
    nix = 1000.0,
    ui0 = 300.0,
    uix = 100.0
)

@variables x t
Dt = Differential(t)
Dx = Differential(x)

nn_manufactured = MMS_CONSTS.nn0 + MMS_CONSTS.nnx*cos(2 * π * MMS_CONSTS.n_waves * x / MMS_CONSTS.L)
function nn_manufactured_f(x, MMS_CONSTS)
    MMS_CONSTS.nn0 + MMS_CONSTS.nnx*cos(2 * π * MMS_CONSTS.n_waves * x / MMS_CONSTS.L)
end

ni_manufactured =  MMS_CONSTS.ni0 + MMS_CONSTS.nix*x/MMS_CONSTS.L #MMS_CONSTS.ni0 + MMS_CONSTS.nix*cos(2 * π * MMS_CONSTS.n_waves * x / MMS_CONSTS.L)
function ni_manufactured_f(x, MMS_CONSTS)
    MMS_CONSTS.ni0 + MMS_CONSTS.nix*x/MMS_CONSTS.L # MMS_CONSTS.ni0 + MMS_CONSTS.nix*cos(2 * π * MMS_CONSTS.n_waves * x / MMS_CONSTS.L)
end 

ui_manufactured = MMS_CONSTS.ui0 + MMS_CONSTS.uix*x/MMS_CONSTS.L #2000 - ni_manufactured #MMS_CONSTS.ui0 + MMS_CONSTS.uix*cos(2 * π * MMS_CONSTS.n_waves * x / MMS_CONSTS.L)
function ui_manufactured_f(x, MMS_CONSTS)
    MMS_CONSTS.ui0 + MMS_CONSTS.uix*x/MMS_CONSTS.L #2000 - ni_manufactured_f(x, MMS_CONSTS)#MMS_CONSTS.ui0 + MMS_CONSTS.uix*cos(2 * π * MMS_CONSTS.n_waves * x / MMS_CONSTS.L)
end

RHS_1 = Dt(nn_manufactured) + Dx(nn_manufactured * MMS_CONSTS.un)
RHS_2 = Dt(ni_manufactured) + Dx(ni_manufactured * ui_manufactured)
RHS_3 = Dt(ni_manufactured * ui_manufactured) + Dx(ni_manufactured * ui_manufactured^2 + ni_manufactured*HallThruster.kB*MMS_CONSTS.ion_temperature)

derivs = expand_derivatives.([RHS_1, RHS_2, RHS_3])

RHS_func = build_function(derivs, [x])
mms! = eval(RHS_func[2]) #return [1] as RHS_1 and [2] as RHS_2, mms([3 3])

simulation = (
    ncells = 100,
    propellant = HallThruster.Xenon,
    ncharge = 1,
    MMS = true,
    mms! = mms!, 
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
    saveat = [end_time],
    tspan = (0., end_time),
    dt = 200e-8, #5e-8
    scheme = (
        flux_function = HallThruster.upwind!,
        limiter = identity,
        reconstruct = false
    ),
)

mutable struct Result
    solution #::Vector{Matrix{Float64}}
    z_cells #::Vector{Float64}
    ncells
    u_exa
    errors
    L_inf
    L_1
end

refinements = MMS_CONSTS.refinements
results = Array{Result, 1}(undef, refinements)
simulations_mms = Array{NamedTuple, 1}(undef, refinements)
n_cells = MMS_CONSTS.n_cells_start
_, __, fluid_ranges, ___ = HallThruster.configure_simulation(simulation)

for refinement in 1:refinements
    global simulations_mms[refinement] = merge(simulation, (; ncells = n_cells, dt = MMS_CONSTS.CFL*MMS_CONSTS.L/(MMS_CONSTS.un*n_cells)))
    global n_cells = n_cells*2
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
    results[refinement] = Result(sol.u, z_cells, simulations_mms[refinement].ncells, u_exa, error, [maximum(error[i, :]) for i in 1:size(error)[1]], [Statistics.mean(error[i, :]) for i in 1:size(error)[1]])
    println("N cells: $(simulations_mms[refinement].ncells) ")
    println("CFL number: $(simulations_mms[refinement].dt*MMS_CONSTS.un/((z_cells[end]-z_cells[1])/simulations_mms[refinement].ncells))") #CFL number
end  

function compute_slope(ncells, errors)
    p = Array{Union{Nothing, Float64}}(nothing, length(ncells)-2)
    for i in 1:length(ncells)-2
        p[i] = log(abs(errors[i+2]-errors[i+1])/abs(errors[i+1]-errors[i]))/log(0.5)
    end 
    return Statistics.mean(p)
end

println("L1 error norm continuity neutral $(compute_slope([results[i].ncells for i in 1:length(results)], [results[i].L_1[1] for i in 1:length(results)]))")
println("L1 error norm continuity ion $(compute_slope([results[i].ncells for i in 1:length(results)], [results[i].L_1[2] for i in 1:length(results)]))")
println("L1 error norm momentum ion $(compute_slope([results[i].ncells for i in 1:length(results)], [results[i].L_1[3] for i in 1:length(results)]))")

#println("L_inf error norm $(compute_slope([results[i].ncells for i in 1:length(results)], [results[i].L_inf[1] for i in 1:length(results)]))")

#plot(log.([length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_1[1] for i in 1:length(results)]))

