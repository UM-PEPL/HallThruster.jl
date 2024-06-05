function allocate_arrays(grid, fluids, anom_model = HallThruster.NoAnom())
    # Number of variables in the state vector U
    nvariables = 0
    for fluid in fluids
        if fluid.type == _ContinuityOnly
            nvariables += 1
        elseif fluid.type == _IsothermalEuler
            nvariables += 2
        elseif fluid.type == _EulerEquations
            nvariables += 3
        end
    end

    ncells = grid.ncells + 2
    nedges = grid.ncells + 1

    U = zeros(nvariables, ncells)
    B = zeros(ncells)
    Aϵ = Tridiagonal(ones(ncells-1), ones(ncells), ones(ncells-1)) #for energy
    bϵ = zeros(ncells) #for energy
    νan = zeros(ncells)
    νc = zeros(ncells)
    νei = zeros(ncells)
    νen = zeros(ncells)
    radial_loss_frequency = zeros(ncells)
    νew_momentum = zeros(ncells)
    νiw = zeros(ncells)
    κ = zeros(ncells)
    μ = zeros(ncells)
    ϕ = zeros(ncells)
    ∇ϕ = zeros(ncells)
    ne = zeros(ncells)
    nϵ = zeros(ncells)
    Tev = zeros(ncells)
    pe = zeros(ncells)
    ∇pe = zeros(ncells)
    ue = zeros(ncells)
    F = zeros(nvariables, nedges)
    UL = zeros(nvariables, nedges)
    UR = zeros(nvariables, nedges)
    Z_eff = zeros(ncells)
    λ_global = zeros(length(fluids))
    νe = zeros(ncells)
    νiz = zeros(ncells)
    νex = zeros(ncells)
    K   = zeros(ncells)
    ji  = zeros(ncells)

    ohmic_heating = zeros(ncells)
    wall_losses = zeros(ncells)
    inelastic_losses = zeros(ncells)

    ncharge = maximum(f.species.Z for f in fluids)
    ni = zeros(ncharge, ncells)
    ui = zeros(ncharge, ncells)
    niui = zeros(ncharge, ncells)
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
    dt_iz = zeros(ncells)
    dt_E = zeros(ncells)
    dt_u = zeros(nedges)
    dt_cell = zeros(ncells)
    dt = zeros(1)

    cache = (;
                Aϵ, bϵ, nϵ, B, νan, νc, μ, ϕ, ∇ϕ, ne, Tev, pe, ue, ∇pe,
                νen, νei, radial_loss_frequency, νew_momentum, νiw, νe, κ, F, UL, UR, Z_eff, λ_global, νiz, νex, K, Id, ji,
                ni, ui, Vs, niui, nn, k, u1, γ_SEE, cell_cache_1,
                error_integral, Id_smoothed, anom_multiplier, smoothing_time_constant,
                errors, dcoeffs,
                ohmic_heating, wall_losses, inelastic_losses,
                channel_area, dA_dz, channel_height, inner_radius, outer_radius, tanδ,
                anom_variables, dt_iz, dt_E, dt_u, dt, dt_cell
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
"""
function run_simulation(
        config::Config;
        grid::HallThrusterGrid = EvenGrid(0),
        ncells = 0,
        dt, duration, nsave,  restart = nothing,
        CFL = 0.799, adaptive = false,
        control_current = false, target_current = 0.0,
        Kp = 0.0, Ti = Inf, Td = 0.0, time_constant = 5e-4,
        dtmin = 0.0, dtmax = Inf,
        verbose = true,
    )

    #check that Landmark uses the correct thermal conductivity
    if config.LANDMARK && !(config.conductivity_model isa LANDMARK_conductivity)
        error("LANDMARK configuration needs to use the LANDMARK thermal conductivity model.")
    end

    # If duration and/or dt are provided with units, convert to seconds and then strip units
    duration = convert_to_float64(duration, u"s")
    dt = convert_to_float64(dt, u"s")

    fluids, fluid_ranges, species, species_range_dict, is_velocity_index = configure_fluids(config)

    if (grid.ncells == 0)
        # backwards compatability for unspecified grid type
        grid1d = generate_grid(config.thruster.geometry, config.domain, EvenGrid(ncells))
    else
        grid1d = generate_grid(config.thruster.geometry, config.domain, grid)
    end

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
        U, cache = load_restart(grid1d, fluids, config, restart)
    else
        U, cache = allocate_arrays(grid1d, fluids, config.anom_model)
    end

    tspan = (0.0, duration)
    saveat = range(tspan[1], tspan[2], length = nsave)

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
    cache.dt_cell .= dt

    # Simulation parameters
    params = (;
        ncells = ncells+2,
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
        dtmin, dtmax,
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

    # Set up
    sol_info = @timed solve(U, params, tspan; saveat)

    sol = sol_info.value
    sim_time = sol_info.time

    # Print some diagnostic information
    if sol.retcode != :Success && verbose
        println("Simulation exited at t = $(sol.t[end]) with retcode :$(sol.retcode) in $(sim_time) seconds.")
    end

    return sol
end
