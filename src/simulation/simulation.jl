function allocate_arrays(grid, config)
    (;ncharge, anom_model) = config

    # made less general to handle common use cases as part of fluid refactor
    nvariables = 1 + 2 * ncharge    # 1 variable for neutrals and 2 for each ion fluid

    ncells = grid.ncells + 2
    nedges = grid.ncells + 1

    # Main state vector
    U = zeros(nvariables, ncells)

    # Caches for energy solve
    Aϵ = Tridiagonal(ones(ncells-1), ones(ncells), ones(ncells-1))
    bϵ = zeros(ncells)

    # Collision frequencies
    νan = zeros(ncells)
    νc = zeros(ncells)
    νei = zeros(ncells)
    νen = zeros(ncells)
    radial_loss_frequency = zeros(ncells)
    νew_momentum = zeros(ncells)
    νiw = zeros(ncells)
    νe = zeros(ncells)
    νiz = zeros(ncells)
    νex = zeros(ncells)

    # Magnetic field
    B = zeros(ncells)

    # Conductivity and mobility
    κ = zeros(ncells)
    μ = zeros(ncells)

    # Potential and electric field
    ϕ = zeros(ncells)
    ∇ϕ = zeros(ncells)

    # Electron number density
    ne = zeros(ncells)

    # Electron energy density
    nϵ = zeros(ncells)

    # Electron temperature and energy [eV]
    Tev = zeros(ncells)
    ϵ = zeros(ncells)

    # Electron pressure and pressure gradient
    pe = zeros(ncells)
    ∇pe = zeros(ncells)

    # Electron axial velocity and kinetic energy
    ue = zeros(ncells)
    K = zeros(ncells)

    λ_global = zeros(ncharge + 1)

    # Electron source terms
    ohmic_heating = zeros(ncells)
    wall_losses = zeros(ncells)
    inelastic_losses = zeros(ncells)

    # Effective charge number
    Z_eff = zeros(ncells)

    # Ion density, velocity, and number flux
    ni = zeros(ncharge, ncells)
    ui = zeros(ncharge, ncells)
    niui = zeros(ncharge, ncells)

    # ion current
    ji  = zeros(ncells)

    # Neutral density
    nn = zeros(ncells)
    γ_SEE = zeros(ncells)
    Id  = [0.0]
    error_integral = [0.0]
    Id_smoothed = [0.0]
    Vs = [0.0]
    anom_multiplier = [1.0]
    smoothing_time_constant = [0.0]
    errors = [0.0, 0.0, 0.0]
    dcoeffs = [0.0, 0.0, 0.0, 0.0]

    # Edge state caches
    F = zeros(nvariables, nedges)
    UL = zeros(nvariables, nedges)
    UR = zeros(nvariables, nedges)

    # timestepping caches
    k = copy(U)
    u1 = copy(U)

    # other caches
    cell_cache_1 = zeros(ncells)

    # Plume divergence variables
    channel_area = zeros(ncells)    # Area of channel / plume
    dA_dz = zeros(ncells)           # derivative of area w.r.t. axial coordinate
    channel_height = zeros(ncells)  # Height of channel / plume (outer - inner)
    inner_radius = zeros(ncells)    # Channel/plume inner radius
    outer_radius = zeros(ncells)    # Channel/plume outer radius
    tanδ = zeros(ncells)            # Tangent of divergence half-angle

    # Anomalous transport variables
    anom_variables = allocate_anom_variables(anom_model, size(U, 2))

    # Timesteps
    dt_iz = zeros(1)
    dt = zeros(1)
    dt_E = zeros(1)
    dt_u = zeros(nedges)

    cache = (;
        Aϵ, bϵ, nϵ, B, νan, νc, μ, ϕ, ∇ϕ, ne, ϵ, Tev, pe, ue, ∇pe,
        νen, νei, radial_loss_frequency, νew_momentum, νiw, νe, κ, F, UL, UR, Z_eff, λ_global, νiz, νex, K, Id, ji,
        ni, ui, Vs, niui, nn, k, u1, γ_SEE, cell_cache_1,
        error_integral, Id_smoothed, anom_multiplier, smoothing_time_constant,
        errors, dcoeffs,
        ohmic_heating, wall_losses, inelastic_losses,
        channel_area, dA_dz, channel_height, inner_radius, outer_radius, tanδ,
        anom_variables, dt_iz, dt_E, dt_u, dt
    )

    return U, cache
end

function setup_simulation(
        config::Config;
        grid::HallThrusterGrid = EvenGrid(0),
        ncells = 0,
        dt, restart = nothing,
        CFL = 0.799, adaptive = false,
        control_current = false, target_current = 0.0,
        Kp = 0.0, Ti = Inf, Td = 0.0, time_constant = 5e-4,
        dtmin = 1e-10, dtmax = 1e-7, max_small_steps = 100,
    )

    #check that Landmark uses the correct thermal conductivity
    if config.LANDMARK && !(config.conductivity_model isa LANDMARK_conductivity)
        error("LANDMARK configuration needs to use the LANDMARK thermal conductivity model.")
    end

    # If dt is provided with units, convert to seconds and then strip units
    dt = convert_to_float64(dt, u"s")
    dtbase = dt

    fluids, fluid_ranges, species, species_range_dict, is_velocity_index = configure_fluids(config)

    if (grid.ncells == 0)
        # backwards compatability for unspecified grid type
        grid1d = generate_grid(config.thruster.geometry, config.domain, EvenGrid(ncells))
    else
        grid1d = generate_grid(config.thruster.geometry, config.domain, grid)
    end

    # load collisions and reactions
    directories = config.reaction_rate_directories
    if config.LANDMARK
        landmark_rates_dir = joinpath(LANDMARK_FOLDER, "reactions")
        directories = [landmark_rates_dir; directories]
    end
    ionization_reactions = load_ionization_reactions(config.ionization_model, species; directories)
    ionization_reactant_indices = reactant_indices(ionization_reactions, species_range_dict)
    ionization_product_indices = product_indices(ionization_reactions, species_range_dict)

    excitation_reactions = load_excitation_reactions(config.excitation_model, species; directories)
    excitation_reactant_indices = reactant_indices(excitation_reactions, species_range_dict)

    electron_neutral_collisions = load_elastic_collisions(config.electron_neutral_model, species; directories)

    index = configure_index(fluids, fluid_ranges)

    use_restart = restart !== nothing

    if use_restart
        U, cache = load_restart(grid1d, config, restart)
    else
        U, cache = allocate_arrays(grid1d, config)
    end

    z_cell = grid1d.cell_centers
    z_edge = grid1d.edges

    # Fill up cell lengths and magnetic field vectors
    for (i, z) in enumerate(grid1d.cell_centers)
        cache.B[i] = config.thruster.magnetic_field(z)
    end

    Δz_cell, Δz_edge = grid_spacing(grid1d)

    mi = config.propellant.m

    # make the adaptive timestep independent of input condition
    if adaptive
        dt = 100 * eps() # small initial timestep to initialize everything

        # force the CFL to be no higher than 0.799 for adaptive timestepping
        # this limit is mainly due to empirical testing, but there
        # may be an analytical reason the ionization timestep cannot use a CFL >= 0.8
        if CFL >= 0.8
            @warn("CFL for adaptive timestepping set higher than stability limit of 0.8. Setting CFL to 0.799.")
            CFL = 0.799
        end
    end

    cache.smoothing_time_constant[] = time_constant
    cache.dt .= dt

    # Simulation parameters
    params = (;
        ncells = grid1d.ncells+2,
        ncharge = config.ncharge,
        mi,
        config = config,
        ϕ_L = config.discharge_voltage + config.cathode_potential,
        ϕ_R = config.cathode_potential,
        Te_L = config.anode_Te,
        Te_R = config.cathode_Te,
        L_ch = config.thruster.geometry.channel_length,
        A_ch = config.thruster.geometry.channel_area,
        z_cell,
        z_edge,
        index, cache, fluids, fluid_ranges, species_range_dict, is_velocity_index,
        iteration = [-1],
        ionization_reactions,
        ionization_reactant_indices,
        ionization_product_indices,
        excitation_reactions,
        excitation_reactant_indices,
        electron_neutral_collisions,
        dt = [dt],
        CFL,
        adaptive,
        background_neutral_velocity = background_neutral_velocity(config),
        background_neutral_density = background_neutral_density(config),
        Bmax = maximum(cache.B),
        γ_SEE_max = 1 - 8.3 * sqrt(me / mi),
        Δz_cell, Δz_edge,
        control_current, target_current, Kp, Ti, Td,
        exit_plane_index = findfirst(>=(config.thruster.geometry.channel_length), z_cell) - 1,
        dtbase, dtmin, dtmax, max_small_steps,
        # landmark benchmark uses pe = 3/2 ne Te, otherwise use pe = ne Te
        pe_factor = config.LANDMARK ? 3/2 : 1.0
    )

    # Compute maximum allowed iterations
    if !use_restart
        initialize!(U, params)

        # Initialize the anomalous collision frequency using a two-zone Bohm approximation for the first iteration
        TwoZoneBohm(1//160, 1//16)(params.cache.νan, params)

        # Initialize plume
        update_plume_geometry!(U, params, initialize = true)
    end

    # make values in params available for first timestep
    update_heavy_species!(U, params)
    update_electrons!(params)

    return U, params
end

"""
    $(SIGNATURES)
Given a state vector `U` and params struct generated by `setup_simulation`, run for provided `duration`, saving `nsave` snapshots.
"""
function run_from_setup(U, params; duration, nsave, verbose = true)
    duration = convert_to_float64(duration, u"s")

    tspan = (0.0, duration)
    saveat = range(tspan[1], tspan[2], length = nsave)

    sol_info = @timed solve(U, params, tspan; saveat)
    sol = sol_info.value

    # Print some diagnostic information
    if sol.retcode != :Success && verbose
        println("Simulation exited at t = $(sol.t[end]) with retcode :$(sol.retcode) in $(sol_info.time) seconds.")
    end

    return sol
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
"""
function run_simulation(config::Config; duration, nsave, verbose = true, kwargs...)
    U, params = setup_simulation(config; kwargs...)
    return run_from_setup(U, params; duration, nsave, verbose)
end
