struct HyperbolicScheme{F,L}
    flux_function::F  # in-place flux function
    limiter::L # limiter
    reconstruct::Bool
end

Base.@kwdef mutable struct MultiFluidSimulation{IC,B1,B2,B3,B4,S,F,L,CB,SP,BP} #could add callback, or autoselect callback when in MMS mode
    grid::Grid1D
    fluids::Vector{Fluid}     # An array of user-defined fluids.
    # This will give us the capacity to more easily do shock tubes (and other problems)
    # without Hall thruster baggage
    initial_condition::IC
    boundary_conditions::Tuple{B1,B2,B3,B4}   # Tuple of left and right boundary conditions, subject to the approval of PR #10
    end_time::Float64    # How long to simulate
    scheme::HyperbolicScheme{F,L} # Flux, Limiter
    source_term!::S  # Source term function. This can include reactons, electric field, and MMS terms
    source_potential!::SP #potential source term
    boundary_potential!::BP #boundary conditions potential
    saveat::Vector{Float64} #when to save
    timestepcontrol::Tuple{Float64,Bool} #sets timestep (first argument) if second argument false. if second argument (adaptive) true, given dt is ignored.
    callback::CB
    solve_energy::Bool
    OVS::Verification
end

function update_heavy_species!(dU, U, params, t) #get source and BCs for potential from params
    ####################################################################
    #extract some useful stuff from params
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    index = params.index

    F, UL, UR, Q = params.cache.F, params.cache.UL, params.cache.UR, params.cache.Q
    B = params.cache.B

    z_cell, z_edge, cell_volume = params.z_cell, params.z_edge, params.cell_volume
    scheme = params.scheme

    ncells = size(U, 2) - 2

    mi = params.mi
    un = fluids[1].conservation_laws.u

    #=
    #fluid BCs
    @views apply_bc!(U[1:index.lf, :], params.BCs[1], :left, params.Te_L, mi)
    @views apply_bc!(U[1:index.lf, :], params.BCs[2], :right, params.Te_R, mi)
    =#

    # Inject anode mass flow at left boundary
    U[index.ρn, begin] = params.mdot_a / params.A_ch / un

    # Neumann BC for neutral density at right boundary
    U[index.ρn, end] = U[index.ρn, end-1]

    # Ion boundary conditions
    for i in 1:params.ncharge
        # Additional neutral flux at left boundary due to recombined ions
        U[index.ρn, 1] -= U[index.ρiui[i], begin] / un

        u_bohm = sqrt(2/3 * i * e * U[index.Tev] / mi)

        # Make sure ions at left boundary satisfy Bohm condition
        boundary_flux = U[index.ρiui[i], begin+1]
        boundary_velocity = min(-u_bohm, boundary_flux / U[index.ρi[i], begin+1])
        boundary_density = boundary_flux / boundary_velocity
        U[index.ρi[i], begin] = boundary_density
        U[index.ρiui[i], begin] = U[index.ρiui[i], begin+1]

        # Neumann BC for ions at the right boundary
        U[index.ρi[i], end] = U[index.ρi[i], end-1]
        U[index.ρiui[i], end] = U[index.ρiui[i], end-1]
    end

    #fluid computations, electron in implicit
    @views compute_edge_states!(UL[1:index.lf, :], UR[1:index.lf, :], U[1:index.lf, :], scheme)
    @views compute_fluxes!(F[1:index.lf, :], UL[1:index.lf, :], UR[1:index.lf, :], fluids, fluid_ranges, scheme, 0.0.*U[index.pe, :])

    # Compute heavy species source terms
    @inbounds  for i in 2:(ncells + 1)
        @turbo Q .= 0.0

        # source term - includes acceleration, ionization, and other user-speficied source terms
        params.ion_source_term!(Q, U, params, i)

        #Compute dU/dt
        left = left_edge(i)
        right = right_edge(i)

        Δz = z_edge[right] - z_edge[left]

        # Ion and neutral fluxes
        @turbo @views @. dU[1:index.lf, i] = (F[1:index.lf, left] - F[1:index.lf, right]) / Δz + Q[1:index.lf]
        dU[index.nϵ, i] = Q[index.nϵ]

        # Ion diffusion term
        η = params.δ * sqrt(2*e*U[index.Tev, i]/(3*mi))
        @views dU[2:index.lf, i] += η*(U[2:index.lf, i-1] - 2U[2:index.lf, i] + U[2:index.lf, i+1])/(Δz)^2

    end

    if !params.implicit_energy
        update_electron_energy!(dU, U, params, t)
    end

end

function update_electron_energy!(dU, U, params, t)
    #########################################################
    #ELECTRON SOLVE

    index = params.index
    μ = params.cache.μ
    z_cell = params.z_cell
    ncells = size(U, 2) - 2
    F = params.cache.F
    z_edge = params.z_edge

    #=
    apply_bc_electron!(U, params.BCs[3], :left, index)
    apply_bc_electron!(U, params.BCs[4], :right, index)
    =#

    # Dirichlet BCs for the electrons
    U[index.nϵ, 1] = U[index.ne, 1] * params.Te_L
    U[index.nϵ, end] = U[index.ne, end] * params.Te_R

    @inbounds for i in 2:ncells+1

        left = left_edge(i)
        right = right_edge(i)

        # Upwinded first derivatives
        ue = U[index.ue, i]
        ne = U[index.ne, i]

        Δz = z_edge[right] - z_edge[left] #boundary cells same size with upwind
        #ue⁺ = max(ue, 0.0) / ue
        #ue⁻ = min(ue, 0.0) / ue
        ue⁺ = smooth_if(ue, 0.0, 0.0, ue, 0.01) / ue
        ue⁻ = smooth_if(ue, 0.0, ue, 0.0, 0.01) / ue


        advection_term = (
            ue⁺ * (ue * U[index.nϵ, i] - U[index.ue, i-1] * U[index.nϵ, i-1]) -
            ue⁻ * (ue * U[index.nϵ, i] - U[index.ue, i+1] * U[index.nϵ, i+1])
        ) / Δz

        diffusion_term_1 = -(
            ue⁺ * (-(μ[i-1] * U[index.nϵ, i-1] - μ[i] * U[index.nϵ, i]) * (U[index.Tev, i-1] - U[index.Tev, i])) -
            ue⁻ * (μ[i+1] * U[index.nϵ, i+1] - μ[i] * U[index.nϵ, i]) * (U[index.Tev, i+1] - U[index.Tev, i])
        ) / Δz^2

        # Biased differences to account for different cell sizes
        if i == 2
            diffusion_term_2 = μ[i] * U[index.nϵ, i] * (8U[index.Tev, i-1] - 12U[index.Tev, i] + 4U[index.Tev, i+1])
            diffusion_term_2 /= 3*Δz^2
        elseif i == ncells + 1
            diffusion_term_2 = μ[i] * U[index.nϵ, i] * (4U[index.Tev, i-1] - 12U[index.Tev, i] + 8U[index.Tev, i+1])
            diffusion_term_2 /= 3*Δz^2
        else
            # Central second derivatives
            diffusion_term_2 = μ[i] * U[index.nϵ, i] * (U[index.Tev, i-1] - 2U[index.Tev, i] + U[index.Tev, i+1])
            diffusion_term_2 /= Δz^2
        end

        dU[index.nϵ, i] += -5/3 * advection_term + 10/9 * (diffusion_term_1 + diffusion_term_2)
    end

    return nothing
end


function update_values!(U, params)

    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    index = params.index

    B = params.cache.B

    z_cell, z_edge, cell_volume = params.z_cell, params.z_edge, params.cell_volume
    ncells = size(U, 2) - 2

    #update useful quantities relevant for potential, electron energy and fluid solve
    L_ch = params.L_ch
    mi = m(fluids[1])

    @inbounds @views for i in 1:(ncells + 2)
        @views U[index.ne, i] = max(1e13, electron_density(U[:, i], fluid_ranges) / mi)
        U[index.Tev, i] = max(0.1, U[index.nϵ, i]/U[index.ne, i])
        U[index.pe, i] = U[index.nϵ, i]
        params.cache.νan[i] = params.anom_model(i, U, params)
        params.cache.νc[i] = electron_collision_freq(U[index.Tev, i], U[1, i]/mi , U[index.ne, i], mi)
        params.cache.μ[i] = electron_mobility(params.cache.νan[i], params.cache.νc[i], B[i])
        #params.cache.μ[i] = 10.0
    end

    # update electrostatic potential
    solve_potential_edge!(U, params)

    @inbounds @views for i in 1:(ncells + 2)
        @views U[index.grad_ϕ, i] = first_deriv_central_diff_pot(U[index.ϕ, :], params.z_cell, i)
        U[index.ue, i] = electron_velocity(U, params, i)
        #U[index.ue, i] = -100.0
    end

    U[index.ue, 1] = U[index.ue, 2]
    U[index.ue, end] = U[index.ue, end-1]

end

condition(u,t,integrator) = t < 1
function affect!(integrator)
    update_values!(integrator.u, integrator.p)
end

#=
function run_simulation(sim, restart_file = nothing) #put source and Bcs potential in params
    species, fluids, fluid_ranges, species_range_dict = configure_simulation(sim)

    lf = fluid_ranges[end][end]
    index = (;lf = lf, nϵ = lf+1, Tev = lf+2, ne = lf+3, pe = lf+4, ϕ = lf+5, grad_ϕ = lf+6, ue = lf+7)

    use_restart = restart_file !== nothing

    if use_restart
        U, grid, B = read_restart(restart_file)
        _, cache = allocate_arrays(grid, fluids)
        cache.B .= B
    else
        grid = sim.grid
        U, cache = allocate_arrays(grid, fluids)
        initial_condition!(@views(U[1:index.nϵ, :]), grid.cell_centers, sim.initial_condition, fluids)
        precompute_bfield!(cache.B, grid.cell_centers)
    end

    scheme = sim.scheme
    source_term! = sim.source_term!
    timestep = sim.timestepcontrol[1]
    adaptive = sim.timestepcontrol[2]
    tspan = (0.0, sim.end_time)
    anom_model = TwoZoneBohm(1/160, 1/16)

    reactions = load_ionization_reactions(species)
    landmark = load_landmark()
    #ϕ_hallis, grad_ϕ_hallis = load_hallis_for_input()

    BCs = sim.boundary_conditions

    OVS = Verification(0, 0, 0)

    ϕ_L = 300.0
    ϕ_R = 0.0
    L_ch = 0.025
    params = (; L_ch, ϕ_L, ϕ_R, OVS, index, cache, fluids, fluid_ranges, species_range_dict, z_cell=grid.cell_centers,
              z_edge=grid.edges, cell_volume=grid.cell_volume, source_term!, reactions,
              scheme, BCs, dt=timestep, source_potential! = sim.source_potential!,
              boundary_potential! = sim.boundary_potential!, landmark, implicit_energy = sim.solve_energy,
              anom_model,
    )

    #PREPROCESS
    #make values in params available for first timestep
    update_values!(U, params)

    if sim.callback !== nothing
        cb = CallbackSet(DiscreteCallback(condition, affect!, save_positions=(false,false)), sim.callback)
    else
        cb = DiscreteCallback(condition, affect!, save_positions=(false,false))
    end

    implicit_energy = sim.solve_energy

    maxiters = Int(ceil(1000 * tspan[2] / timestep))

    if implicit_energy
        splitprob = SplitODEProblem{true}(update_electron_energy!, update_heavy_species!, U, tspan, params)
        sol = solve(splitprob, KenCarp47(); saveat=sim.saveat, callback=cb,
        adaptive=adaptive, dt=timestep, dtmax = timestep)
    else
        #alg = SSPRK22()
        alg = AutoTsit5(Rosenbrock23())
        prob = ODEProblem{true}(update_heavy_species!, U, tspan, params)
        sol = solve(prob, alg; saveat=sim.saveat, callback=cb,
        adaptive=adaptive, dt=timestep, dtmax=timestep, maxiters = maxiters)
    end

    return HallThrusterSolution(sol, params)
end
=#

function run_simulation(config)
    U, params = configure_simulation(config)
    update_values!(U, params)

    tspan = (0.0, config.simulation_time)
    saveat = LinRange(tspan[1], tspan[2], config.nsave)
    maxiters = Int(ceil(1000 * tspan[2] / config.dtmax))

    cb = nothing
    if cb !== nothing
        callback = CallbackSet(DiscreteCallback(condition, affect!, save_positions=(false,false)), sim.callback)
    else
        callback = DiscreteCallback(condition, affect!, save_positions=(false,false))
    end

    if config.implicit_energy
        splitprob = SplitODEProblem{true}(update_electron_energy!, update_heavy_species!, U, tspan, params)
        sol = solve(splitprob, KenCarp47(); saveat, callback,
        adaptive=config.adaptive, dt=config.dt, dtmax = config.dtmax)
    else
        alg = AutoTsit5(Rosenbrock23())
        prob = ODEProblem{true}(update_heavy_species!, U, tspan, params)
        sol = solve(
            prob, alg; saveat, callback,
            adaptive=config.adaptive,
            dt=config.dt,
            dtmax=config.dtmax,
            maxiters = maxiters
        )
    end

    return HallThrusterSolution(sol, params)
end


function initial_condition!(U, z_cell, IC!, fluids)
    for (i, z) in enumerate(z_cell)
        @views IC!(U[:, i], z, fluids, z_cell[end])
    end
end
