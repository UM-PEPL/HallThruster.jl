@kwdef mutable struct SimParams{G <: GridSpec, C <: CurrentController}
    grid::G
    dt::Float64
    duration::Float64
    # Optional parameters
    num_save::Int        = 1000
    restart::Bool        = false
    CFL::Float64         = 0.799
    adaptive::Bool       = false
    min_dt::Float64      = 1e-10
    max_dt::Float64      = 1e-7
    max_small_steps::Int = 100
    verbose::Bool        = true
    print_errors::Bool   = true
    current_control::C   = NoController()
end

function setup_simulation(config::Config, sim::SimParams;
        postprocess::Postprocess = Postprocess(), restart = "",)
    #check that Landmark uses the correct thermal conductivity
    if config.LANDMARK && !(config.conductivity_model isa LANDMARK_conductivity)
        error("LANDMARK configuration needs to use the LANDMARK thermal conductivity model.")
    end

    fluids, fluid_ranges, species, species_range_dict, is_velocity_index = configure_fluids(config)
    index = configure_index(fluids, fluid_ranges)

    # load collisions and reactions
    ionization_reactions = load_ionization_reactions(
        config.ionization_model, unique(species);
        directories = config.reaction_rate_directories,)
    ionization_reactant_indices = reactant_indices(ionization_reactions, species_range_dict)
    ionization_product_indices = product_indices(ionization_reactions, species_range_dict)

    excitation_reactions = load_excitation_reactions(
        config.excitation_model, unique(species),)
    excitation_reactant_indices = reactant_indices(excitation_reactions, species_range_dict)

    electron_neutral_collisions = load_elastic_collisions(
        config.electron_neutral_model, unique(species),)

    # Generate grid and allocate state
    grid = generate_grid(config.thruster.geometry, config.domain, sim.grid)

    if sim.restart
        U, cache = load_restart(grid, config, restart)
    else
        U, cache = allocate_arrays(grid, config)
    end

    # Fill up cell lengths and magnetic field vectors
    thruster = config.thruster
    itp = LinearInterpolation(thruster.magnetic_field.z, thruster.magnetic_field.B)
    @. cache.B = itp(grid.cell_centers)

    # make the adaptive timestep independent of input condition
    dt = sim.dt
    if sim.adaptive
        dt = 100 * eps() # small initial timestep to initialize everything

        # force the CFL to be no higher than 0.799 for adaptive timestepping
        # this limit is mainly due to empirical testing, but there
        # may be an analytical reason the ionization timestep cannot use a CFL >= 0.8
        if sim.CFL >= 0.8
            if sim.print_errors
                @warn("CFL for adaptive timestepping set higher than stability limit of 0.8. Setting CFL to 0.799.")
            end
            sim.CFL = 0.799
        end
    end

    cache.dt .= dt

    # Simulation parameters
    params = (;
        iteration = [-1],
        dt = [dt],
        grid,
        config = config,
        simulation = sim,
        postprocess,
        # fluid bookkeeping
        index, cache, fluids, fluid_ranges, species_range_dict, is_velocity_index,
        # reactions
        ionization_reactions,
        ionization_reactant_indices,
        ionization_product_indices,
        excitation_reactions,
        excitation_reactant_indices,
        electron_neutral_collisions,
        # Physics stuff
        background_neutral_velocity = background_neutral_velocity(config),
        background_neutral_density = background_neutral_density(config),
        γ_SEE_max = 1 - 8.3 * sqrt(me / config.propellant.m),
        min_Te = min(config.anode_Tev, config.cathode_Tev),
    )

    # Compute maximum allowed iterations
    if !sim.restart
        initialize!(U, params)
        initialize_plume_geometry(params)

        # Initialize the anomalous collision frequency using a two-zone Bohm approximation for the first iteration
        TwoZoneBohm(1 // 160, 1 // 16)(params.cache.νan, params)
    end

    # make values in params available for first timestep
    update_heavy_species!(U, params)
    update_electrons!(params)

    return U, params
end

function setup_simulation(
        config::Config;
        duration, nsave,
        grid::GridSpec = EvenGrid(0),
        ncells = 0,
        dt, restart = nothing,
        CFL = 0.799, adaptive = false,
        control_current = false, target_current = 0.0,
        Kp = 0.0, Ti = Inf, Td = 0.0, time_constant = 5e-4,
        dtmin = 1e-10, dtmax = 1e-7, max_small_steps = 100,
        verbose = true, print_errors = true,
)
    # Set up PID controller, if requested
    if control_current
        current_control = PIDController(
            target_value = target_current,
            proportional_constant = Kp,
            integral_constant = Kp / Ti,
            derivative_constant = Kp * Td,
            smoothing_frequency = 1 / time_constant,
        )
    else
        current_control = NoController()
    end

    if (grid.num_cells == 0)
        # backwards compatability for unspecified grid type
        grid = EvenGrid(ncells)
    end

    simulation = SimParams(;
        duration,
        grid,
        dt,
        num_save = nsave,
        CFL,
        restart = restart !== nothing,
        adaptive,
        min_dt = dtmin,
        max_dt = dtmax,
        max_small_steps,
        current_control,
        verbose,
        print_errors,
    )

    return setup_simulation(config, simulation; restart)
end

"""
    $(SIGNATURES)
Given a state vector `U` and params struct generated by `setup_simulation`, run for provided `duration`, saving `nsave` snapshots.
"""
function run_from_setup(U, params)
    tspan = (0.0, params.simulation.duration)
    saveat = range(tspan[1], tspan[2], length = params.simulation.num_save)

    sol_info = @timed solve(U, params, tspan; saveat)
    sol = sol_info.value

    # Print some diagnostic information
    if sol.retcode != :Success && params.simulation.verbose
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
function run_simulation(config::Config; kwargs...)
    U, params = setup_simulation(config; kwargs...)
    return run_from_setup(U, params)
end

function run_simulation(config::Config, sim::SimParams; kwargs...)
    U, params = setup_simulation(config, sim; kwargs...)
    return run_from_setup(U, params)
end
