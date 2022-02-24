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
    mutable struct EnergyOVS
Enables setting mu, ue, Tev and ne to certain values to very electron energy equation
"""
mutable struct EnergyOVS{F1, F2}
    active ::Int64
    μ::Union{Float64, Nothing}
    ue::Union{Float64, Nothing}
    Tev::F1
    ne::F2
end

"""
    mutable struct Verification
is passed to params to identify if OVS is active.
"""

mutable struct Verification{F1, F2}
    potential ::Int64
    fluid ::Int64
    energy ::EnergyOVS{F1, F2}
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

    left_state_elec = 0.0
    right_state_elec = left_state_elec
    BCs_elec = (HallThruster.Dirichlet_energy(left_state_elec), HallThruster.Dirichlet_energy(right_state_elec))

    n_save = 1
    saveat = if n_save == 1
        [MMS_CONSTS.max_end_time]
    else
        LinRange(0.0, MMS_CONSTS.max_end_time, n_save) |> collect
    end

    simulation = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(HallThruster.SPT_100, MMS_CONSTS.n_cells_start), 
    boundary_conditions = (BCs[1], BCs[2], BCs_elec[1], BCs_elec[2]),
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
        #set up U and params
        U = Array{Union{Nothing, Float64}}(nothing, 10, n_cells+2)
        U .= 0.0
        grid = HallThruster.generate_grid(HallThruster.SPT_100, n_cells)
        lf = 3
        index = (;lf = lf, nϵ = lf+1, Tev = lf+2, ne = lf+3, pe = lf+4, ϕ = lf+5, grad_ϕ = lf+6, ue = lf+7)
        U[index.ne, :] .= 1.0 #electron density 
        fluids = [HallThruster.Fluid(HallThruster.Species(HallThruster.Xenon, 0), HallThruster.ContinuityOnly(300.0, 300.0))
        HallThruster.Fluid(HallThruster.Species(HallThruster.Xenon, 1), HallThruster.IsothermalEuler(0.0))]
        L_ch = 0.025
        ϕ_L = 100.0
        ϕ_R = 0.0
        A = Tridiagonal(ones(n_cells), ones(n_cells+1), ones(n_cells))
        b = zeros(n_cells + 1)
        μ = ones(n_cells + 2)
        cache = (; A, b, μ)
        OVS = Verification(1, 0, EnergyOVS(0, nothing, nothing, nothing, nothing))
        params = (; L_ch, ϕ_L, ϕ_R, OVS, index, cache, fluids, z_cell=grid.cell_centers, z_edge = grid.edges)
        HallThruster.solve_potential_edge!(U, params)
        z_cells = params.z_cell
        z_edge = params.z_edge
        ϕ = U[params.index.ϕ, :]
        u_exa = Array{Union{Nothing, Float64}}(nothing, length(ϕ)-3)
        for i in 1:length(ϕ)-3
            u_exa[i] = 25000.0*z_edge[i+1]^2 - 3250.0*z_edge[i+1] + 100.0
        end
        error = abs.(u_exa - ϕ[2:end-2])
        results[refinement] = Result(ϕ, 1.0, z_edge, n_cells-1, u_exa, error, maximum(error), mean(error))
        n_cells = n_cells*2
    end
    return results
end

function compute_slope(refinements, errors)
    p = Array{Union{Nothing, Float64}}(nothing, refinements-2)
    for i in 1:refinements-2
        p[i] = log(abs(errors[i+2]-errors[i+1])/abs(errors[i+1]-errors[i]))/log(0.5)
    end 
    return Statistics.mean(p)
end

function evaluate_slope(results, MMS_CONSTS)
    nfluids = size(results[1].u_exa)[1]
    if nfluids > MMS_CONSTS.refinements
        nfluids = 1
    end
    L_1 = Array{Union{Nothing, Float64}}(nothing, nfluids)
    L_inf = Array{Union{Nothing, Float64}}(nothing, nfluids)
    for j in 1:nfluids
        L_1[j] = compute_slope(MMS_CONSTS.refinements, [results[i].L_1[j] for i in 1:length(results)])
        L_inf[j] = compute_slope(MMS_CONSTS.refinements, [results[i].L_inf[j] for i in 1:length(results)])
    end
    return L_1, L_inf
end


function perform_OVS_elecenergy(; MMS_CONSTS, fluxfn, reconstruct)

    #create a template definition of source! function somewhere
    function source!(Q, U, params, i)
        mms!(@views(Q[1:4]), [params.z_cell[i]])
    end
    
    function IC!(U, z, fluids, L)
        ρ2 = MMS_CONSTS.n0 + MMS_CONSTS.nx
        u1 = MMS_CONSTS.u0
        ρ1 = MMS_CONSTS.n0 + MMS_CONSTS.nx
        Tev = MMS_CONSTS.Tev0 #+ MMS_CONSTS.Tev_elec_max*sin(π * z / (MMS_CONSTS.L))
        ne = (MMS_CONSTS.n0 + MMS_CONSTS.nx) / fluids[1].species.element.m
        U .= SA[ρ1, ρ2, ρ2*u1, ne*Tev]
        return U
    end

    u1 = MMS_CONSTS.u0
    ρ1 = MMS_CONSTS.n0 + MMS_CONSTS.nx
    T1 = MMS_CONSTS.T_constant

    left_state = [ρ1, ρ1, ρ1*u1] # [ρ1, ρ1*u1, ρ1*E] 
    right_state = [ρ1, ρ1, ρ1*(u1+0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Dirichlet(right_state))

    left_state_elec = ρ1/HallThruster.Xenon.m*MMS_CONSTS.Tev0
    right_state_elec = left_state_elec
    BCs_elec = (HallThruster.Dirichlet_energy(left_state_elec), HallThruster.Dirichlet_energy(right_state_elec))
    
    Tev_func(z) = MMS_CONSTS_ELEC.Tev0 #+ MMS_CONSTS_ELEC.Tev_elec_max*sin(2 * π * z / (MMS_CONSTS_ELEC.L))
    ne_func(z) = MMS_CONSTS_ELEC.n0 + MMS_CONSTS_ELEC.nx*cos(2 * π * MMS_CONSTS_ELEC.n_waves * z / MMS_CONSTS_ELEC.L)

    verification = Verification(0, 1, EnergyOVS(1, MMS_CONSTS.μ,  MMS_CONSTS.ue, Tev_func, ne_func))

    n_save = 100
    saveat = if n_save == 1
        [MMS_CONSTS.max_end_time]
    else
        LinRange(0.0, MMS_CONSTS.max_end_time, n_save) |> collect
    end

    cb_convergence = DiffEqCallbacks.TerminateSteadyState(1e-8, 1e-6, DiffEqCallbacks.allDerivPass)
    cb = cb_convergence

    simulation = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(HallThruster.SPT_100, MMS_CONSTS.n_cells_start), 
    boundary_conditions = (BCs[1], BCs[2], BCs_elec[1], BCs_elec[2]),
    scheme = HallThruster.HyperbolicScheme(fluxfn, HallThruster.minmod, reconstruct),
    initial_condition = IC!, 
    source_term! = source!,
    source_potential! = nothing,
    boundary_potential! = nothing, 
    fluids = [HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.ContinuityOnly(MMS_CONSTS.u_constant, MMS_CONSTS.T_constant));
    HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 1), HallThruster.IsothermalEuler(MMS_CONSTS.T_constant))],
    #[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())], 
    end_time = MMS_CONSTS.max_end_time, 
    saveat = saveat,
    timestepcontrol = (1e-6, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
    callback = cb, 
    solve_energy = false,
    verification = verification
    )

    #insert config
    config = (
        anode_potential = 300.0,
        cathode_potential = 0.0,
        anode_Te = 3.0,
        cathode_Te = 3.0,
        restart_file = nothing,
        radial_loss_coefficients = (0.4, 1.0),
        wall_collision_frequencies = (1e7, 0.0),
        geometry = HallThruster.SPT_100,
        anode_mass_flow_rate = 5e-6,
        neutral_velocity = MMS_CONSTS.u_constant,
        neutral_temperature = MMS_CONSTS.T_constant,
        ion_diffusion_coeff = 0.0,
        implicit_energy = false,
        propellant = HallThruster.Xenon,
        ncharge = 1,
        verification = verification,
        solve_ion_energy = false,
        ion_temperature = 1000.0,
        anom_model = HallThruster.TwoZoneBohm(1/160, 1/16),
        energy_equation = :LANDMARK,
        ionization_coeffs = :LANDMARK,
        electron_pressure_coupled = false,
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
        sol = HallThruster.run_simulation(simulation, config)
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
