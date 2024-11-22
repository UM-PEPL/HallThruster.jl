function setup_sim_state(config::Config, sim::SimParams)
    config = validate_config(config)

    fluids, fluid_ranges, species, species_range_dict, is_velocity_index = configure_fluids(config)
    index = configure_index(fluids, fluid_ranges)

    # Load collisions and reactions
    directories = config.reaction_rate_directories
    if config.LANDMARK
        landmark_rates_dir = joinpath(LANDMARK_FOLDER, "reactions")
        pushfirst!(directories, landmark_rates_dir)
    end
    ionization_reactions = load_ionization_reactions(
        config.ionization_model, species; directories,)
    ionization_reactant_indices = reactant_indices(ionization_reactions, species_range_dict)
    ionization_product_indices = product_indices(ionization_reactions, species_range_dict)

    excitation_reactions = load_excitation_reactions(
        config.excitation_model, species; directories,)
    excitation_reactant_indices = reactant_indices(excitation_reactions, species_range_dict)

    electron_neutral_collisions = load_elastic_collisions(
        config.electron_neutral_model, species; directories,)

    # Generate grid
    thruster = config.thruster
    geom = thruster.geometry
    grid = generate_grid(geom, config.domain, sim.grid)
    U, cache = allocate_arrays(grid, config)

    # Evaluate magnetic field at cell centers
    bfield = thruster.magnetic_field
    bfield_func = LinearInterpolation(bfield.z, bfield.B)
    for (i, z) in enumerate(grid.cell_centers)
        cache.B[i] = bfield_func(z)
    end

    # make the adaptive timestep independent of input condition
    dt = sim.dt
    cache.dt .= dt

    if sim.adaptive
        dt = 100 * eps() # small initial timestep to initialize everything

        # force the CFL to be no higher than 0.799 for adaptive timestepping
        # this limit is mainly due to empirical testing, but there
        # may be an analytical reason the ionization timestep cannot use a CFL >= 0.8
        if sim.CFL >= 0.8
            @warn("CFL for adaptive timestepping set higher than stability limit of 0.8. Setting CFL to 0.799.")
            sim.CFL = 0.799
        end
    end

    exit_plane_index = findfirst(>=(geom.channel_length), grid.cell_centers) - 1

    state = (;
        config,
        # redundant with stuff in geometry
        L_ch = thruster.geometry.channel_length,
        A_ch = thruster.geometry.channel_area,
        exit_plane_index,
        # just use the grid
        ncells = grid.num_cells,
        z_cell = grid.cell_centers,
        z_edge = grid.edges,
        Δz_cell = grid.dz_cell,
        Δz_edge = grid.dz_edge,
        # physical properties (unsure if we need them)
        background_neutral_velocity = background_neutral_velocity(config),
        background_neutral_density = background_neutral_density(config),
        γ_SEE_max = 1 - 8.3 * sqrt(me / config.propellant.m),
        Te_min = min(config.anode_Tev, config.cathode_Tev),
        # fluid stuff
        index,
        cache,
        fluids,
        fluid_ranges,
        species_range_dict,
        is_velocity_index,
        iteration = [-1],
        # reactions
        ionization_reactions,
        ionization_reactant_indices,
        ionization_product_indices,
        excitation_reactions,
        excitation_reactant_indices,
        electron_neutral_collisions,
        # just use simparams
        dt = [dt],
        CFL = sim.CFL,
        adaptive = sim.adaptive,
        dt_base = sim.dt,
        min_dt = sim.min_dt,
        max_dt = sim.max_dt,
        max_small_steps = sim.max_small_steps,
        current_control = sim.current_control,
    )

    # Compute maximum allowed iterations
    initialize!(U, state)

    # Initialize the anomalous collision frequency using a two-zone Bohm approximation for the first iteration
    TwoZoneBohm(1 // 160, 1 // 16)(state.cache.nu_anom, state)

    # Initialize plume
    update_plume_geometry!(U, state; initialize = true)
    #end

    # make values in params available for first timestep
    update_heavy_species!(U, state)
    update_electrons!(state)

    return U, state
end
