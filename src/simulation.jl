function allocate_arrays(grid, fluids) #rewrite allocate arrays as function of set of equations, either 1, 2 or 3
    # Number of variables in the state vector U
    nvariables = 0
    for i in 1:length(fluids)
        if fluids[i].conservation_laws.type == _ContinuityOnly
            nvariables += 1
        elseif fluids[i].conservation_laws.type == _IsothermalEuler
            nvariables += 2
        elseif fluids[i].conservation_laws.type == _EulerEquations
            nvariables += 3
        end
    end

    ncells = grid.ncells + 2
    nedges = grid.ncells + 1

    U = zeros(nvariables + 1, ncells)
    A = Tridiagonal(ones(nedges-1), ones(nedges), ones(nedges-1)) #for potential
    b = zeros(nedges) #for potential equation
    B = zeros(ncells)
    Aϵ = Tridiagonal(ones(ncells-1), ones(ncells), ones(ncells-1)) #for energy
    bϵ = zeros(ncells) #for energy
    νan = zeros(ncells)
    νc = zeros(ncells)
    νei = zeros(ncells)
    νen = zeros(ncells)
    νw = zeros(ncells)
    μ = zeros(ncells)
    ϕ = zeros(nedges)
    ϕ_cell = zeros(ncells)
    ∇ϕ = zeros(ncells)
    ne = zeros(ncells)
    Tev = zeros(ncells)
    pe = zeros(ncells)
    ∇pe = zeros(ncells)
    ue = zeros(ncells)
    F = zeros(nvariables+1, nedges)
    UL = zeros(nvariables+1, nedges)
    UR = zeros(nvariables+1, nedges)

    cache = (; A, b, Aϵ, bϵ, B, νan, νc, μ, ϕ, ϕ_cell, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen, νei, νw, F, UL, UR)
    return U, cache
end

function initial_condition!(U, z_cell, IC!, fluids)
    for (i, z) in enumerate(z_cell)
        @views IC!(U[:, i], z, fluids, z_cell[end])
    end
end

function run_simulation(config, timestep, end_time, ncells, nsave; alg = SSPRK22(;stage_limiter!, step_limiter! = stage_limiter!), restart_file = nothing)

    fluids, fluid_ranges, species, species_range_dict = configure_fluids(config)

    index = configure_index(fluid_ranges)

    use_restart = restart_file !== nothing

    if use_restart
        U, grid, B = read_restart(config.restart_file)
        _, cache = allocate_arrays(grid, fluids)
        cache.B .= B
    else
        grid = generate_grid(config.thruster.geometry, ncells, config.domain)
        U, cache = allocate_arrays(grid, fluids)
        initial_condition!(@views(U[1:index.nϵ, :]), grid.cell_centers, config.initial_condition!, fluids)
    end

    tspan = (0.0, end_time)
    saveat = LinRange(tspan[1], tspan[2], nsave)

    ionization_reactions = load_ionization_reactions(config.ionization_model, species)

    for (i, z) in enumerate(grid.cell_centers)
        cache.B[i] = config.thruster.magnetic_field(z)
    end

    # callback for calling the update_values! function at each timestep
    update_callback = DiscreteCallback(Returns(true), update_values!, save_positions=(false,false))

    # Choose which cache variables to save and set up saving callback
    fields_to_save = (:μ, :Tev, :ϕ, :∇ϕ, :ne, :pe, :ue, :∇pe, :νan, :νc, :νen, :νei, :νw, :ϕ_cell)

    function save_func(u, t, integrator)
        (; μ, Tev, ϕ, ∇ϕ, ne, pe, ue, ∇pe, νan, νc, νen, νei, νw, ϕ_cell) = integrator.p.cache
        return deepcopy((; μ, Tev, ϕ, ∇ϕ, ne, pe, ue, ∇pe, νan, νc, νen, νei, νw, ϕ_cell))
    end

    saved_values = SavedValues(
        Float64, NamedTuple{fields_to_save, NTuple{length(fields_to_save), Vector{Float64}}}
    )
    saving_callback = SavingCallback(save_func, saved_values; saveat)

    niters = round(Int, tspan[2] / timestep)

    progress_bar = make_progress_bar(niters, timestep, config)

    # Assemble callback set
    if config.callback !== nothing
        callbacks = CallbackSet(update_callback, saving_callback, config.callback)
    else
        callbacks = CallbackSet(update_callback, saving_callback)
    end

    # Simulation parameters
    params = (;
        config = config,
        ϕ_L = config.discharge_voltage + config.cathode_potential,
        ϕ_R = config.cathode_potential,
        Te_L = config.anode_Te,
        Te_R = config.cathode_Te,
        L_ch = config.thruster.geometry.channel_length,
        A_ch = config.thruster.geometry.channel_area,
        reactions = ionization_reactions,
        z_cell=grid.cell_centers,
        z_edge=grid.edges,
        dt=timestep,
        progress_bar,
        index, cache, fluids, fluid_ranges, species_range_dict,
        iteration = [-1]
    )

    # Compute maximum allowed iterations
    maxiters = Int(ceil(1000 * tspan[2] / timestep))

    #make values in params available for first timestep
    update_values!(U, params)

    # Set up ODE problem and solve
    prob = ODEProblem{true}(update_heavy_species!, U, tspan, params)
	sol = try
        solve(
            prob, alg; saveat, callback=callbacks,
            adaptive=false, dt=timestep, dtmax=10*timestep, dtmin = timestep/10, maxiters = maxiters,
	    )
    catch e
        stop_progress_bar!(progress_bar, params)
        println("There was an error")
        throw(e)
    end

    # Print some diagnostic information
    if sol.retcode == :NaNDetected
        println("Simulation failed with NaN detected at t = $(sol.t[end])")
    elseif sol.retcode == :InfDetected
        println("Simulation failed with Inf detected at t = $(sol.t[end])")
    end

    stop_progress_bar!(progress_bar, params)

    # Return the solution
    return HallThrusterSolution(sol, params, saved_values.saveval)
end