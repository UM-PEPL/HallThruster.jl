@public run_simulation

function setup_simulation(
        config::Config, sim::SimParams;
        postprocess::Union{Postprocess, Nothing} = nothing, include_dirs = String[],
        restart::AbstractString = "",
    )

    # check that Landmark uses the correct thermal conductivity
    if config.LANDMARK && !(config.conductivity_model isa LANDMARK_conductivity)
        error("LANDMARK configuration needs to use the LANDMARK thermal conductivity model.")
    end

    # We arrange the fluid containers in three different structures for convenience.
    # First, an array of (;continuity, isothermal) for each propellant species.
    # This is used in initialization and in computing boundary conditions.
    fluids_by_propellant = [allocate_fluids(propellant, sim.grid.num_cells) for propellant in config.propellants]

    # Second, a single NamedTuple of (;continuity, isothermal) for all propellants.
    # This is used in the convective update.
    continuity = vcat([x.continuity for x in fluids_by_propellant]...)
    isothermal = vcat([x.isothermal for x in fluids_by_propellant]...)
    fluid_containers = (; continuity, isothermal)

    # Finally, a single flat array of fluid containers, which we use for reaction calculations.
    fluid_array = vcat([[fluid.continuity..., fluid.isothermal...] for fluid in fluids_by_propellant]...)
    species = unique([fl.species for fl in fluid_array])


    # Load reactions and collisions either from file or generate from species list
    ei_reactions, excitation_reactions, electron_neutral_collisions = load_reactions(
        config.propellant_config, species,
        config.ionization_model, config.excitation_model, config.electron_neutral_model;
        directories = config.reaction_rate_directories
    )

    # Get reactant and product indices
    ei_reactant_indices = reactant_indices(ei_reactions, fluid_array)
    ei_product_indices = product_indices(ei_reactions, fluid_array)
    excitation_reactant_indices = reactant_indices(excitation_reactions, fluid_array)
    electron_neutral_indices = reactant_indices(electron_neutral_collisions, fluid_array)

    # Generate grid and allocate state
    grid = generate_grid(sim.grid, config.thruster.geometry, config.domain)
    cache = allocate_arrays(grid, config)

    # Load magnetic field
    thruster = config.thruster
    load_magnetic_field!(thruster.magnetic_field; include_dirs)

    # Fill up cell lengths and magnetic field vectors
    itp = LinearInterpolation(thruster.magnetic_field.z, thruster.magnetic_field.B)
    @. cache.B = itp(grid.cell_centers) * config.magnetic_field_scale

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
    # IMPORTANT: for effective precompilation, we do not include any types
    # in here that are non-concrete or likely to change between runs.
    # For instance, any method that accepts a Config{...} will be recompiled for
    # all definitions of {...}, regardless of whether it actually needs to be
    # Therefore, we do not include `config` in the `params` struct, instead
    # passing it around sparingly to the methods that really need it.
    # Except for `sim`, nothing in this struct should have type parameters.
    # For convenience, the method `params_from_config` copies concretely-typed
    # values from `config` and reinserts them into params.
    params = (;
        # non-concretely-typed, changes based on run, requires recompilation
        # TODO: is this true anymore?
        params_from_config(config)...,
        # concretely-typed except for PID controller, not too bad
        simulation = sim,
        # Remainder is concretely-typed
        iteration = [-1],
        dt = [dt],
        grid,
        postprocess = if isnothing(postprocess)
            Postprocess()
        else
            postprocess
        end,
        # fluid bookkeeping - concretely-typed
        cache,
        # reactions - concretely-typed
        ei_reactions,
        ei_reactant_indices,
        ei_product_indices,
        excitation_reactions,
        excitation_reactant_indices,
        electron_neutral_collisions,
        electron_neutral_indices,
        min_Te = min(config.anode_Tev, config.cathode_Tev),
        fluid_containers,
        fluid_array,
        fluids_by_propellant,
    )

    # Initialize ion and electron variables
    initialize!(params, config)

    # Initialize the anomalous collision frequency using a
    # two-zone Bohm approximation for the first iteration
    TwoZoneBohm(1 / 160, 1 / 16)(params.cache.νan, params, config)

    if !isempty(restart)
        # Initialize the solution from a restart JSON file
        initialize_from_restart!(params, restart)
    end

    # make values in params available for first timestep
    initialize_plume_geometry(params)
    update_heavy_species!(params)
    update_electrons!(params, config)

    return params
end

function setup_simulation(
        config::Config;
        duration, nsave,
        grid::GridSpec = EvenGrid(0),
        ncells = 0,
        dt, restart::AbstractString = "",
        CFL = 0.799, adaptive = false,
        control_current = false, target_current = 0.0,
        Kp = 0.0, Ti = Inf, Td = 0.0, time_constant = 5.0e-4,
        dtmin = 1.0e-10, dtmax = 1.0e-7, max_small_steps = 100,
        verbose = true, print_errors = true,
    )
    dt = convert_to_float64(dt, units(:s))
    duration = convert_to_float64(duration, units(:s))
    dtmin = convert_to_float64(dtmin, units(:s))
    dtmax = convert_to_float64(dtmax, units(:s))

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
Given a params struct generated by `setup_simulation`, run for provided `duration`, saving `nsave` snapshots.
"""
function run_from_setup(params, config)
    tspan = (0.0, params.simulation.duration)
    saveat = range(tspan[1], tspan[2], length = params.simulation.num_save)

    sol_info = @timed solve(params, config, tspan; saveat)
    sol = sol_info.value

    # Print some diagnostic information
    if sol.retcode != :Success && params.simulation.verbose
        @info "Simulation exited at t = $(sol.t[end]) with retcode :$(sol.retcode) in $(sol_info.time) seconds."
    end

    return sol
end

"""
	$(TYPEDSIGNATURES)
Run a Hall thruster simulation using the provided `Config` and `SimParams` objects.
Returns a `Solution` object.

## Arguments
- `config::Config`: contains geometry, plasma properties, and numerical information about the simulation. See [Configuration](@ref) for more information.
- `sim::SimParams`: contains grid generation and timestepping information. See [Simulations](@ref) for more information.
- `postprocess::Union{Postprocess, Nothing}`: contains file to which output is to be written and specifies what kind of output to write. If `nothing`, no output is written to file. See [Postprocessing](@ref) for more information.
- `restart`::String:An optional path to a JSON file containing plasma data. If non-empty, the solution will be restarted from that file.
"""
function run_simulation(config::Config, sim::SimParams; postprocess = nothing, restart::String = "", kwargs...)
    params = setup_simulation(config, sim; postprocess, restart, kwargs...)
    return run_from_setup(params, config)
end

"""
    $(TYPEDSIGNATURES)
**Deprecated**. Please use `run_simulation(::Config, ::SimParams; kwargs...)`

Run a Hall thruster simulation using the provided Config object.

## Arguments
- `config`: a `Config` containing simulation parameters.
- `dt`: The timestep, in seconds. Typical values are O(10 ns) (1e-8 seconds).
- `duration`: How long to run the simulation, in seconds (simulation time, not wall time). Typical runtimes are O(1 ms) (1e-3 seconds).
- `ncells`: How many cells to use. Typical values are 100 - 1000 cells.
- `nsave`: How many frames to save.
"""
function run_simulation(config::Config; kwargs...)
    params = setup_simulation(config; kwargs...)
    return run_from_setup(params, config)
end
