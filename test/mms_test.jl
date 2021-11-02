using Test, Documenter, HallThruster, StaticArrays, BenchmarkTools, Symbolics, DifferentialEquations, Statistics, Plots

configure_simulation = HallThruster.configure_simulation
generate_grid = HallThruster.generate_grid
inlet_neutral_density = HallThruster.inlet_neutral_density
load_ionization_reactions = HallThruster.load_ionization_reactions
reconstruct! = HallThruster.reconstruct!
compute_fluxes! = HallThruster.compute_fluxes!
left_edge = HallThruster.left_edge
right_edge = HallThruster.right_edge


Te_func = z -> 30 * exp(-(2(z - SPT_100.channel_length) / 0.033)^2)
ϕ_func = z -> 300 * (1 - 1/(1 + exp(-1000 * (z - SPT_100.channel_length))))
ni_func = z -> 2000 #1e6

end_time = 1e-5 #30e-5

#need to merge MMS_CONSTS and simulation parameters to not have a discrepancy between the two
#also for ions calculate the velocity, this is not fixed and not given, ie wrong for now
#next check convergence

simulation = (
    ncells = 100,
    propellant = HallThruster.Xenon,
    ncharge = 1,
    geometry = HallThruster.SPT_100,
    neutral_temperature = 500.,
    neutral_velocity = 300.,
    ion_temperature = 0.0,
    initial_Te = Te_func,
    initial_ϕ = ϕ_func,
    initial_ni = ni_func,
    solve_Te = false,
    solve_ne = false,
    inlet_mdot = 5e-6,
    saveat = [end_time],
    tspan = (0., end_time),
    dt = 5e-8, #5e-8
    scheme = (
        flux_function = HallThruster.upwind!,
        limiter = identity,
        reconstruct = false
    ),
)

const MMS_CONSTS = (
    n_cells_start = 10,
    refinements = 2,
    n_waves = 2.0,
    un = 300, 
    L = 0.05,
    ion_temperature = 500,
    nn0 = 1000.0,
    nnx = 1000.0,
    ni0 = 2000.0,
    nix = 0000.0,
    ui0 = 300.0,
    uix = 0.0
)

@variables x t
Dt = Differential(t)
Dx = Differential(x)

nn_manufactured = MMS_CONSTS.nn0 + MMS_CONSTS.nnx*cos(2 * π * MMS_CONSTS.n_waves * x / MMS_CONSTS.L)
function nn_manufactured_f(x, MMS_CONSTS)
    MMS_CONSTS.nn0 + MMS_CONSTS.nnx*cos(2 * π * MMS_CONSTS.n_waves * x / MMS_CONSTS.L)
end

ni_manufactured =  MMS_CONSTS.ni0 + MMS_CONSTS.nix*cos(2 * π * MMS_CONSTS.n_waves * x / MMS_CONSTS.L)
function ni_manufactured_f(x, MMS_CONSTS)
    MMS_CONSTS.ni0 + MMS_CONSTS.nix*cos(2 * π * MMS_CONSTS.n_waves * x / MMS_CONSTS.L)
end 

ui_manufactured = 2000 - ni_manufactured #MMS_CONSTS.ui0 + MMS_CONSTS.uix*cos(2 * π * MMS_CONSTS.n_waves * x / MMS_CONSTS.L)
function ui_manufactured_f(x, MMS_CONSTS)
    2000 - ni_manufactured_f(x, MMS_CONSTS)#MMS_CONSTS.ui0 + MMS_CONSTS.uix*cos(2 * π * MMS_CONSTS.n_waves * x / MMS_CONSTS.L)
end

RHS_1 = Dt(nn_manufactured) + Dx(nn_manufactured * MMS_CONSTS.un)
RHS_2 = Dt(ni_manufactured) + Dx(ni_manufactured * ui_manufactured)
RHS_3 = Dt(ni_manufactured * ui_manufactured) + Dx(ni_manufactured * ui_manufactured^2 + ni_manufactured*HallThruster.kB*simulation.ion_temperature)

derivs = expand_derivatives.([RHS_1, RHS_2, RHS_3])

RHS_func = build_function(derivs, [x])
mms = eval(RHS_func[1]) #return [1] as RHS_1 and [2] as RHS_2, mms([3 3])

#############################################################################################

function update_MMS!(dU, U, params, t) #added mms terms to dU and got rid of source term stuff
	fluids, fluid_ranges = params.fluids, params.fluid_ranges

	F, UL, UR, Q = params.cache

	z_cell, z_edge = params.z_cell, params.z_edge
	scheme = params.scheme

	nvariables = size(U, 1)
    ncells = size(U, 2) - 2

	Te_index = nvariables-2
	ne_index = nvariables-1
	ϕ_index = nvariables

    # zero variables and apply boundary conditions
    # TEMPORARY: apply Neumann BC on right edge
    for i in 1:nvariables
	    dU[i, 1] = 0.0
	    dU[i, ncells+1] = 0.0
        U[i, end] = U[i, end-1]
        #U[i, end] = U[i, 1]
        #U[i, 1] = U[i, end]
    end

    reconstruct!(UL, UR, U, scheme)
    compute_fluxes!(F, UL, UR, fluids, fluid_ranges, scheme)

    # deleted all the source term parts and electric field stuff

    for i in 2:ncells+1 #since more faces than cells
    # Compute dU/dt
		left = left_edge(i) #i-1
		right = right_edge(i) #i

		Δz = z_edge[right] - z_edge[left]

        @views dU[1, i] = (F[1, left] - F[1, right])/Δz + params.mms([z_cell[i] z_cell[i]])[1] #added mms terms here

		for (i, j) in enumerate(fluid_ranges[2:end])
			@views dU[j[1], i] = (F[j[1], left] - F[j[1], right])/Δz #+ params.mms([z_cell[i] z_cell[i]])[j[1]] #added mms terms here
            @views dU[j[2], i] = (F[j[2], left] - F[j[2], right])/Δz #+ params.mms([z_cell[i] z_cell[i]])[j[2]] #added mms terms here
		end
    end
    return nothing
end


#should be fine, just need to make it run multiple times and record the values,

function run_simulation_MMS(sim, mms) #added mms as input and to params

    species, fluids, fluid_ranges, species_range_dict = configure_simulation(sim)
    z_cell, z_edge = generate_grid(sim.geometry, sim.ncells)

    U, cache = allocate_arrays_MMS(sim)

    initial_condition_MMS!(U, z_cell, sim, fluid_ranges)

    scheme = sim.scheme

    #reactions = load_ionization_reactions(species)

    params = (; # deleted reactions
        cache,
        fluids,
        fluid_ranges,
        species_range_dict,
        z_cell,
        z_edge,
        scheme,
        mms
    )

    prob = ODEProblem{true}(update_MMS!, U, sim.tspan, params)
    sol = solve(prob, Tsit5(), dt = sim.dt, saveat = sim.saveat) #SSPRK33
    return sol
end

function allocate_arrays_MMS(sim) # got rid of electrons and electric field variables
    # Number of variables in the state vector U
    # U = [nn, ni1, ni1ui1..., niN, niNuiN]
    nvariables = 1 + 2 * sim.ncharge

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

function initial_condition_MMS!(U, z_cell, sim, fluid_ranges) #got rid of electron and electric field variables
    nvariables = size(U, 1)
    nn = inlet_neutral_density(sim)
    un = sim.neutral_velocity

    nn_index = 1

    for (i, z) in enumerate(z_cell)
        if i < length(z_cell)/2
            U[nn_index, i] = 2000.0
        else
            U[nn_index, i] = 0.0
        end
        
        ni = sim.initial_ni(z)
        # ions initialized with equal densities, same velocity as neutrals
        for j in fluid_ranges[2:end] #make same init condition
            n_index = j[1]
            nu_index = j[2]
            if i < length(z_cell)/2
                U[n_index, i] = ni*sin(z)
                U[nu_index, i] = ni * un * sin(z)
            else 
                U[n_index, i] = ni*sin(z)
                U[nu_index, i] = ni * un*sin(z)
            end
            
        end
    end

    return U
end

##############################################################################################################################################

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
_, __, fluid_ranges, ___ = configure_simulation(simulation)

for refinement in 1:refinements
    global simulations_mms[refinement] = merge(simulation, (; ncells = n_cells))
    global n_cells = n_cells*2
end

for refinement in 1:refinements
    sol = run_simulation_MMS(simulations_mms[refinement], mms)
    z_cells = HallThruster.generate_grid(simulation.geometry, simulations_mms[refinement].ncells)[1]
    u_exa = Array{Union{Nothing, Float64}}(nothing, length(sol.u[1][:, 1]), length(z_cells))
    for (i, z_cell) in enumerate(z_cells)
        u_exa[1, i] = nn_manufactured_f(z_cell, MMS_CONSTS)
        for (j, index) in enumerate(fluid_ranges[2:end])
            u_exa[index[1], i] = ni_manufactured_f(z_cell, MMS_CONSTS)
            u_exa[index[2], i] = ni_manufactured_f(z_cell, MMS_CONSTS) * ui_manufactured_f(z_cell, MMS_CONSTS)
        end
    end
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

println("L1 error norm $(compute_slope([results[i].ncells for i in 1:length(results)], [results[i].L_1[1] for i in 1:length(results)]))")
#println("L_inf error norm $(compute_slope([results[i].ncells for i in 1:length(results)], [results[i].L_inf[1] for i in 1:length(results)]))")


plot(log.([length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_1[1] for i in 1:length(results)]))

