
@kwdef struct SimParams{G <: HallThrusterGrid, C <: CurrentController}
    # Grid setup
    grid::G = EvenGrid(0)
    ncells::Int = 0

    # Timestep control
    base_dt::Float64
    min_dt::Float64 = 1e-8
    max_dt::Float64 = 1e-7
    CFL::Float64 = 0.799
    adaptive::Bool = false
    max_small_steps::Int = 100

    # PID control
    current_control::C = NoController()

    # Reporting
    verbose::Bool = true
    show_errors::Bool = true
end

function setup_simulation(
        config::Config;
        grid::HallThrusterGrid = EvenGrid(0),
        ncells = 0,
        dt,
        restart = nothing,
        CFL = 0.799,
        adaptive = false,
        control_current = false,
        target_current = 0.0,
        Kp = 0.0,
        Ti = Inf,
        Td = 0.0,
        time_constant = 5e-4,
        dtmin = 1e-10,
        dtmax = 1e-7,
        max_small_steps = 100,
)
    #check that Landmark uses the correct thermal conductivity
    if config.LANDMARK && !(config.conductivity_model isa LANDMARK_conductivity)
        error(
            "LANDMARK configuration needs to use the LANDMARK thermal conductivity model."
        )
    end

    # If dt is provided with units, convert to seconds and then strip units
    dtbase = dt

    fluids, fluid_ranges, species, species_range_dict, is_velocity_index = configure_fluids(
        config
    )

    if (grid.ncells == 0)
        # backwards compatability for unspecified grid type
        grid_1d = generate_grid(config.thruster.geometry, config.domain, EvenGrid(ncells))
    else
        grid_1d = generate_grid(config.thruster.geometry, config.domain, grid)
    end

    # load collisions and reactions
    directories = config.reaction_rate_directories
    if config.LANDMARK
        landmark_rates_dir = joinpath(LANDMARK_FOLDER, "reactions")
        directories = [landmark_rates_dir; directories]
    end
    ionization_reactions = load_ionization_reactions(
        config.ionization_model, species; directories,
    )
    ionization_reactant_indices = reactant_indices(ionization_reactions, species_range_dict)
    ionization_product_indices = product_indices(ionization_reactions, species_range_dict)

    excitation_reactions = load_excitation_reactions(
        config.excitation_model, species; directories,
    )
    excitation_reactant_indices = reactant_indices(excitation_reactions, species_range_dict)

    electron_neutral_collisions = load_elastic_collisions(
        config.electron_neutral_model, species; directories,
    )

    index = configure_index(fluids, fluid_ranges)

    U, cache = allocate_arrays(grid_1d, config)

    z_cell = grid_1d.cell_centers
    z_edge = grid_1d.edges

    # Fill up cell lengths and magnetic field vectors
    thruster = config.thruster
    bfield = thruster.magnetic_field
    bfield_func = LinearInterpolation(bfield.z, bfield.B)
    for (i, z) in enumerate(grid_1d.cell_centers)
        cache.B[i] = bfield_func(z)
    end

    Δz_cell, Δz_edge = grid_spacing(grid_1d)

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
        ncells = grid_1d.ncells + 2,
        ncharge = config.ncharge,
        mi,
        config = config,
        V_L = config.discharge_voltage,
        V_R = config.cathode_potential,
        Te_L = config.anode_Tev,
        Te_R = config.cathode_Tev,
        Te_min = min(config.anode_Tev, config.cathode_Tev),
        L_ch = thruster.geometry.channel_length,
        A_ch = thruster.geometry.channel_area,
        z_cell,
        z_edge,
        index,
        cache,
        fluids,
        fluid_ranges,
        species_range_dict,
        is_velocity_index,
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
        Δz_cell,
        Δz_edge,
        control_current,
        target_current,
        Kp,
        Ti,
        Td,
        exit_plane_index = findfirst(>=(thruster.geometry.channel_length), z_cell) - 1,
        dtbase,
        dtmin,
        dtmax,
        max_small_steps,
        # landmark benchmark uses pe = 3/2 ne Te, otherwise use pe = ne Te
        pe_factor = config.LANDMARK ? 3 / 2 : 1.0,
    )

    # Compute maximum allowed iterations
    initialize!(U, params)

    # Initialize the anomalous collision frequency using a two-zone Bohm approximation for the first iteration
    TwoZoneBohm(1 // 160, 1 // 16)(params.cache.nu_anom, params)

    # Initialize plume
    update_plume_geometry!(U, params; initialize = true)
    #end

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
    tspan = (0.0, duration)
    saveat = range(tspan[1], tspan[2]; length = nsave)

    sol_info = @timed solve(U, params, tspan; saveat)
    sol = sol_info.value

    # Print some diagnostic information
    if sol.retcode != :Success && verbose
        println(
            "Simulation exited at t = $(sol.t[end]) with retcode :$(sol.retcode) in $(sol_info.time) seconds.",
        )
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
    config = validate_config(config)
    U, params = setup_simulation(config; kwargs...)
    return run_from_setup(U, params; duration, nsave, verbose)
end
