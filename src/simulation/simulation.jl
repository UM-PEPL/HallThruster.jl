function allocate_arrays(grid, fluids) #rewrite allocate arrays as function of set of equations, either 1, 2 or 3
    # Number of variables in the state vector U
    nvariables = 0
    for fluid in fluids
        if fluid.conservation_laws.type == _ContinuityOnly
            nvariables += 1
        elseif fluid.conservation_laws.type == _IsothermalEuler
            nvariables += 2
        elseif fluid.conservation_laws.type == _EulerEquations
            nvariables += 3
        end
    end

    ncells = grid.ncells + 2
    nedges = grid.ncells + 1

    U = zeros(nvariables + 1, ncells)
    B = zeros(ncells)
    Aϵ = Tridiagonal(ones(ncells-1), ones(ncells), ones(ncells-1)) #for energy
    bϵ = zeros(ncells) #for energy
    νan = zeros(ncells)
    νc = zeros(ncells)
    νei = zeros(ncells)
    νen = zeros(ncells)
    νew = zeros(ncells)
    νiw = zeros(ncells)
    μ = zeros(ncells)
    ϕ = zeros(ncells)
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

    num_neutral_fluids = count(f -> f.species.Z == 0, fluids)
    ncharge = maximum(f.species.Z for f in fluids)
    ni = zeros(ncharge, ncells)
    ui = zeros(ncharge, ncells)
    niui = zeros(ncharge, ncells)
    nn = zeros(num_neutral_fluids, ncells)
    nn_tot = zeros(ncells)
    Id  = [0.0]
    Vs = [0.0]

    cache = (;
                Aϵ, bϵ, B, νan, νc, μ, ϕ, ∇ϕ, ne, Tev, pe, ue, ∇pe,
                νen, νei, νew, νiw, νe, F, UL, UR, Z_eff, λ_global, νiz, νex, K, Id, ji,
                ni, ui, Vs, niui, nn, nn_tot
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
    restart = nothing, CFL = 1.0, adaptive = false)

    # If duration and/or dt are provided with units, convert to seconds and then strip units
    duration = convert_to_float64(duration, u"s")
    dt = convert_to_float64(dt, u"s")

    fluids, fluid_ranges, species, species_range_dict = configure_fluids(config)
    num_neutral_fluids = count(f -> f.species.Z == 0, fluids)
    grid = generate_grid(config.thruster.geometry, ncells, config.domain)

    # load collisions and reactions
    ionization_reactions = _load_reactions(config.ionization_model, unique(species))
    ionization_reactant_indices = reactant_indices(ionization_reactions, species_range_dict)
    ionization_product_indices = product_indices(ionization_reactions, species_range_dict)

    excitation_reactions = _load_reactions(config.excitation_model, unique(species))
    excitation_reactant_indices = reactant_indices(excitation_reactions, species_range_dict)

    electron_neutral_collisions = _load_reactions(config.electron_neutral_model, unique(species))

    index = configure_index(fluids, fluid_ranges)

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
    fields_to_save() = (
        :μ, :Tev, :ϕ, :∇ϕ, :ne, :pe, :ue, :∇pe, :νan, :νc, :νen,
        :νei, :νew, :νiz, :νex, :νe, :Id, :ni, :ui, :ji, :niui, :nn, :nn_tot
    )
    extract_from_cache(to_save, cache) = NamedTuple{to_save}(cache)

    function save_func(u, t, integrator)
        return deepcopy(extract_from_cache(fields_to_save(), integrator.p.cache))
    end

    saved_values = SavedValues(
        Float64, typeof(extract_from_cache(fields_to_save(), cache))
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
        max_timestep = [dt],
        CFL,
        adaptive,
        num_neutral_fluids,
        background_neutral_velocity = background_neutral_velocity(config),
        background_neutral_density = background_neutral_density(config),
        Bmax = maximum(cache.B),
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

function run_simulation(json_path::String)
    (;design, simulation, parameters) = JSON3.read(read(json_path, String))

    geometry = Geometry1D(;
        channel_length = design.channel_length,
        outer_radius = design.outer_radius,
        inner_radius = design.inner_radius,
    )

    bfield_data = readdlm(design.magnetic_field_file, ',')
    bfield_func = HallThruster.LinearInterpolation(bfield_data[:, 1], bfield_data[:, 2])

    thruster = HallThruster.Thruster(;
        name = design.thruster_name,
        geometry = geometry,
        magnetic_field = bfield_func,
        shielded = design.magnetically_shielded
    )

    propellant = eval(Symbol(design.propellant))
    wall_material = eval(Symbol(design.wall_material))

    anom_model = if simulation.anom_model == "NoAnom"
        NoAnom()
    elseif simulation.anom_model == "Bohm"
        (x -> Bohm(x[1]))(parameters.anom_model_coeffs)
    elseif simulation.anom_model == "TwoZoneBohm"
        (x -> TwoZoneBohm(x[1], x[2]))(parameters.anom_model_coeffs)
    elseif simulation.anom_model == "MultiLogBohm"
        MultiLogBohm(parameters.anom_model_coeffs)
    end

    flux_function = eval(Symbol(simulation.flux_function))
    limiter = eval(Symbol(simulation.limiter))

    config = HallThruster.Config(;
        thruster,
        propellant,
        anom_model,
        domain = (0.0, simulation.cathode_location_m),
        discharge_voltage = design.anode_potential - design.cathode_potential,
        anode_mass_flow_rate = design.anode_mass_flow_rate,
        cathode_potential = design.cathode_potential,
        ncharge = simulation.max_charge,
        wall_loss_model = WallSheath(wall_material, parameters.sheath_loss_coefficient),
        ion_wall_losses = simulation.ion_wall_losses,
        cathode_Te = parameters.cathode_electron_temp_eV,
        LANDMARK = false,
        ion_temperature = parameters.ion_temp_K,
        neutral_temperature = parameters.neutral_temp_K,
        neutral_velocity = parameters.neutral_velocity_m_s,
        electron_ion_collisions = simulation.electron_ion_collisions,
        min_electron_temperature = parameters.cathode_electron_temp_eV,
        transition_function = LinearTransition(parameters.inner_outer_transition_length_m, 0.0),
        electron_pressure_coupled = false,
        scheme = HyperbolicScheme(;
            flux_function, limiter, reconstruct = simulation.reconstruct, WENO = false
        ),
        solve_background_neutrals           = simulation.solve_background_neutrals,
        background_pressure                 = simulation.background_pressure_Torr,
        background_neutral_temperature      = simulation.background_temperature_K,
    )

    solution = run_simulation(config; ncells = simulation.num_cells,
        nsave = simulation.num_save, duration = simulation.duration_s, dt = simulation.dt_s)

    return solution
end
