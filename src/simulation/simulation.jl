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
    Z_eff = zeros(ncells)
    λ_global = zeros(length(fluids))
    νe = zeros(ncells)
    νiz = zeros(ncells)
    νex = zeros(ncells)
    K   = zeros(ncells)
    ji  = zeros(ncells)
    Id  = [0.0]
    cache = (;
                A, b, Aϵ, bϵ, B, νan, νc, μ, ϕ, ϕ_cell, ∇ϕ, ne, Tev, pe, ue, ∇pe,
                νen, νei, νw, νe, F, UL, UR, Z_eff, λ_global, νiz, νex, K, Id, ji
            )

    return U, cache
end

"""
    $(SIGNATURES)
Run a Hall thruster simulation using the provided Config object.

## Arguments
- `config`: a `Config` containing simulation parameters. 
- `dt`: The timestep, in seconds. Typical values are O(10 ns) (1e-8 seconds).
- `duration`: How long to run the simulation, in seconds (simulation time, not wall time). Typical runtimes are O(1 ms) (1e-3 seconds).
- `ncells`: How many cells to use. Typical values are 100 - 1000 cells.
- `nsave`: How many frames to save.
- `alg`: Explicit timestepping algorithm to use. Defaults to second-order strong-stability preserving Runge-Kutta with a stage limiter and a step limiter(`SSPRK22(stage_limiter = stage_limiter!)`), but will work with any explicit timestepping algorithm provided by `OrdinaryDiffEq`.
- `restart`: path to restart file or a HallThrusterSolution object. Defaults to `nothing`.
"""
function run_simulation(config::Config;
    dt, duration, ncells, nsave, alg = SSPRK22(;stage_limiter!, step_limiter! = stage_limiter!),
    restart = nothing, electron_energy_order = 2, CFL = 1.0, adaptive = false)

    fluids, fluid_ranges, species, species_range_dict = configure_fluids(config)
    grid = generate_grid(config.thruster.geometry, ncells, config.domain)

    # load collisions and reactions
    ionization_reactions = _load_reactions(config.ionization_model, species)
    ionization_reactant_indices = reactant_indices(ionization_reactions, species_range_dict)
    ionization_product_indices = product_indices(ionization_reactions, species_range_dict)

    excitation_reactions = _load_reactions(config.excitation_model, species)
    excitation_reactant_indices = reactant_indices(excitation_reactions, species_range_dict)

    electron_neutral_collisions = _load_reactions(config.electron_neutral_model, species)

    index = configure_index(fluid_ranges)

    use_restart = restart !== nothing

    if use_restart
        U, cache = load_restart(grid, fluids, config, restart)
    else
        U, cache = allocate_arrays(grid, fluids)
    end

    tspan = (0.0, duration)
    saveat = LinRange(tspan[1], tspan[2], nsave)

    for (i, z) in enumerate(grid.cell_centers)
        cache.B[i] = config.thruster.magnetic_field(z)
    end

    # callback for calling the update_values! function at each timestep
    update_callback = DiscreteCallback(Returns(true), update_values!, save_positions=(false,false))

    # Choose which cache variables to save and set up saving callback
    fields_to_save = (:μ, :Tev, :ϕ, :∇ϕ, :ne, :pe, :ue, :∇pe, :νan, :νc, :νen, :νei, :νw, :νiz, :νex, :νe, :ϕ_cell, :Id)

    function save_func(u, t, integrator)
        (; μ, Tev, ϕ, ∇ϕ, ne, pe, ue, ∇pe, νan, νc, νen, νei, νw, νiz, νex, νe, ϕ_cell, Id) = integrator.p.cache
        return deepcopy((; μ, Tev, ϕ, ∇ϕ, ne, pe, ue, ∇pe, νan, νc, νen, νei, νw, νiz, νex, νe, ϕ_cell, Id))
    end

    saved_values = SavedValues(
        Float64, NamedTuple{fields_to_save, NTuple{length(fields_to_save), Vector{Float64}}}
    )
    saving_callback = SavingCallback(save_func, saved_values; saveat)

    # Assemble callback set
    if config.callback !== nothing
        callbacks = CallbackSet(update_callback, saving_callback, config.callback)
    else
        callbacks = CallbackSet(update_callback, saving_callback)
    end

    # Simulation parameters
    params = (;
        ncharge = config.ncharge,
        mi = config.propellant.m,
        config = config,
        ϕ_L = config.discharge_voltage + config.cathode_potential,
        ϕ_R = config.cathode_potential,
        Te_L = config.anode_Te,
        Te_R = config.cathode_Te,
        L_ch = config.thruster.geometry.channel_length,
        A_ch = config.thruster.geometry.channel_area,
        z_cell=grid.cell_centers,
        z_edge=grid.edges,
        dt,
        index, cache, fluids, fluid_ranges, species_range_dict,
        iteration = [-1],
        ionization_reactions,
        ionization_reactant_indices,
        ionization_product_indices,
        excitation_reactions,
        excitation_reactant_indices,
        electron_neutral_collisions,
        electron_energy_order,
        max_timestep = [dt],
        CFL,
        adaptive
    )

    # Compute maximum allowed iterations
    maxiters = Int(ceil(1000 * tspan[2] / dt))

    if !use_restart
        initialize!(U, params)
    end

    #make values in params available for first timestep
    update_values!(U, params)

    # Set up ODE problem and solve
    prob = ODEProblem{true}(update_heavy_species!, U, tspan, params)
	sol = solve(
            prob, alg; saveat, callback=callbacks,
            dt=dt, dtmax=10*dt, dtmin = dt/100, maxiters = maxiters,
	    )

    # Print some diagnostic information
    if sol.retcode == :NaNDetected
        println("Simulation failed with NaN detected at t = $(sol.t[end])")
    elseif sol.retcode == :InfDetected
        println("Simulation failed with Inf detected at t = $(sol.t[end])")
    end

    # Return the solution
    return Solution(sol, params, saved_values.saveval)
end