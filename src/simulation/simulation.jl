function allocate_arrays(grid, fluids, anom_model = HallThruster.NoAnom())
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

    ohmic_heating = zeros(ncells)
    wall_losses = zeros(ncells)
    inelastic_losses = zeros(ncells)

    num_neutral_fluids = count(f -> f.species.Z == 0, fluids)
    ncharge = maximum(f.species.Z for f in fluids)
    ni = zeros(ncharge, ncells)
    ui = zeros(ncharge, ncells)
    niui = zeros(ncharge, ncells)
    nn = zeros(num_neutral_fluids, ncells)
    nn_tot = zeros(ncells)
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
                Aϵ, bϵ, B, νan, νc, μ, ϕ, ∇ϕ, ne, Tev, pe, ue, ∇pe,
                νen, νei, νew, νiw, νe, F, UL, UR, Z_eff, λ_global, νiz, νex, K, Id, ji,
                ni, ui, Vs, niui, nn, nn_tot, k, u1, γ_SEE, cell_cache_1,
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

    # If duration and/or dt are provided with units, convert to seconds and then strip units
    duration = convert_to_float64(duration, u"s")
    dt = convert_to_float64(dt, u"s")

    fluids, fluid_ranges, species, species_range_dict = configure_fluids(config)
    num_neutral_fluids = count(f -> f.species.Z == 0, fluids)

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
    saveat = LinRange(tspan[1], tspan[2], nsave)

    z_cell = grid1d.cell_centers
    z_edge = grid1d.edges

    # Fill up cell lengths and magnetic field vectors

    for (i, z) in enumerate(grid1d.cell_centers)
        cache.B[i] = config.thruster.magnetic_field(z)
    end

    Δz_cell, Δz_edge = grid_spacing(grid1d)

    mi = config.propellant.m

    #make the adaptive timestep independent of input condition 
    if adaptive
        dt = 100 * eps()#small initial timestep to initialize everything

        #force the CFL to be no higher than 0.799 for adaptive timestepping
        #this limit is mainly due to empirical testing, but there 
        #may be an analytical reason the ionization timestep cannot use a CFL>=0.8
        if CFL >= 0.8
            @warn("CFL for Adaptive Timestepping Set Higher than Stability Limit of 0.8. Setting CFL to 0.799.")
            CFL = 0.799
        end
    end

    cache.smoothing_time_constant[] = time_constant
    cache.dt .= dt
    cache.dt_cell .= dt


    # Simulation parameters
    params = (;
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
        index, cache, fluids, fluid_ranges, species_range_dict,
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
        num_neutral_fluids,
        background_neutral_velocity = background_neutral_velocity(config),
        background_neutral_density = background_neutral_density(config),
        Bmax = maximum(cache.B),
        γ_SEE_max = 1 - 8.3 * sqrt(me / mi),
        Δz_cell, Δz_edge,
        control_current, target_current, Kp, Ti, Td,
        exit_plane_index = findfirst(>=(config.thruster.geometry.channel_length), z_cell) - 1,
        dtmin, dtmax
    )


    # Compute maximum allowed iterations
    if !use_restart
        initialize!(U, params)
    end

    # Initialize the anomalous collision frequency using a two-zone Bohm approximation for the first iteration
    TwoZoneBohm(1//160, 1//16)(params.cache.νan, params)

    # Initialize plume
    update_plume_geometry!(U, params, initialize = true)

    # make values in params available for first timestep
    update_electrons!(U, params)

    # Set up and solve problem, this form is temporary
    prob = ODEProblem(update_heavy_species!, U, tspan, params)
    sol_info = @timed solve(prob; saveat)

    sol = sol_info.value
    sim_time = sol_info.time

    # Print some diagnostic information
    if sol.retcode != :Success && verbose
        println("Simulation exited at t = $(sol.t[end]) with retcode :$(sol.retcode) in $(sim_time) seconds.")
    end

    return sol
end

function run_simulation(json_content::JSON3.Object; single_section = false, nonstandard_keys = false, verbose = true)

    adaptive = false
    pressure_z0 = NaN
    pressure_dz = NaN
    pressure_pstar = NaN
    pressure_alpha = NaN

    if single_section
        if nonstandard_keys
            (;
                # Design
                thruster_name,
                channel_length, inner_radius, outer_radius,
                magnetic_field_file, magnetically_shielded,
                propellant_material, wall_material,
                anode_potential, cathode_potential,
                anode_mass_flow_rate,
                # Simulation
                anom_model, cathode_location_m, max_charge,
                num_cells, dt_s, duration_s, num_save,
                flux_function, limiter, reconstruct,
                ion_wall_losses, electron_ion_collisions, solve_background_neutrals,
                # Parameters
                anom_model_coeffs, sheath_loss_coefficient,
                ion_temp_K, neutral_temp_K, neutral_velocity_m_s,
                cathode_electron_temp_eV, inner_outer_transition_length_m,
                background_pressure_Torr, background_temperature_K,
                pressure_z0, pressure_dz, pressure_pstar, pressure_alpha
            ) = json_content

            propellant = propellant_material
        else
            (;
                # Design
                thruster_name,
                channel_length, inner_radius, outer_radius,
                magnetic_field_file, magnetically_shielded,
                propellant, wall_material,
                anode_potential, cathode_potential,
                anode_mass_flow_rate,
                # Simulation
                anom_model, cathode_location_m, max_charge,
                num_cells, dt_s, duration_s, num_save,
                flux_function, limiter, reconstruct,
                ion_wall_losses, electron_ion_collisions, solve_background_neutrals,
                # Parameters
                anom_model_coeffs, sheath_loss_coefficient,
                ion_temp_K, neutral_temp_K, neutral_velocity_m_s,
                cathode_electron_temp_eV, inner_outer_transition_length_m,
                background_pressure_Torr, background_temperature_K,
                pressure_z0, pressure_dz, pressure_pstar, pressure_alpha
            ) = json_content
        end

        # Handle optional keys
        if haskey(json_content, :adaptive)
            adaptive = json_content.adaptive
        end

        # Optional parameters for ShiftedTwoZoneBohm
        if anom_model == "ShiftedTwoZone" || anom_model == "ShiftedTwoZoneBohm"
            (;pressure_z0, pressure_dz, pressure_pstar, pressure_alpha) = json_content
        end
    else
        (;design, simulation, parameters) = json_content
        (;
            thruster_name,
            channel_length, inner_radius, outer_radius,
            magnetic_field_file, magnetically_shielded,
            propellant, wall_material,
            anode_potential, cathode_potential,
            anode_mass_flow_rate,
        ) = design

        (;
            anom_model, cathode_location_m, max_charge,
            num_cells, dt_s, duration_s, num_save,
            flux_function, limiter, reconstruct,
            ion_wall_losses, electron_ion_collisions, solve_background_neutrals,
        ) = simulation

        (;
            anom_model_coeffs, sheath_loss_coefficient,
            ion_temp_K, neutral_temp_K, neutral_velocity_m_s,
            cathode_electron_temp_eV, inner_outer_transition_length_m,
            background_pressure_Torr, background_temperature_K,
        ) = parameters

        # Handle optional keys
        if (haskey(simulation, :adaptive))
            adaptive = simulation.adaptive
        end

        # Optional parameters for ShiftedTwoZoneBohm
        if anom_model == "ShiftedTwoZone" || anom_model == "ShiftedTwoZoneBohm"
            (;pressure_z0, pressure_dz, pressure_pstar, pressure_alpha) = parameters
        end
    end

    geometry = Geometry1D(;channel_length, outer_radius, inner_radius)

    bfield_data = readdlm(magnetic_field_file, ',')
    bfield_func = HallThruster.LinearInterpolation(bfield_data[:, 1], bfield_data[:, 2])

    thruster = HallThruster.Thruster(;
        name = thruster_name,
        geometry = geometry,
        magnetic_field = bfield_func,
        shielded = magnetically_shielded
    )

    anom_model = if anom_model == "NoAnom"
        NoAnom()
    elseif anom_model == "Bohm"
        Bohm(anom_model_coeffs[1])
    elseif anom_model == "TwoZoneBohm"
        TwoZoneBohm((anom_model_coeffs[1], anom_model_coeffs[2]))
    elseif anom_model == "MultiLogBohm"
        MultiLogBohm(anom_model_coeffs)
    elseif anom_model == "ShiftedTwoZone" || anom_model == "ShiftedTwoZoneBohm"
        coeff_tuple = (anom_model_coeffs[1], anom_model_coeffs[2])
        ShiftedTwoZoneBohm(coeff_tuple, pressure_z0, pressure_dz, pressure_pstar, pressure_alpha)
    end

    config = HallThruster.Config(;
        thruster,
        propellant = eval(Symbol(propellant)),
        anom_model,
        domain = (0.0, cathode_location_m),
        discharge_voltage = anode_potential - cathode_potential,
        anode_mass_flow_rate = anode_mass_flow_rate,
        cathode_potential = cathode_potential,
        ncharge = max_charge,
        wall_loss_model = WallSheath(eval(Symbol(wall_material)), sheath_loss_coefficient),
        ion_wall_losses = ion_wall_losses,
        cathode_Te = cathode_electron_temp_eV,
        LANDMARK = false,
        ion_temperature = ion_temp_K,
        neutral_temperature = neutral_temp_K,
        neutral_velocity = neutral_velocity_m_s,
        electron_ion_collisions = electron_ion_collisions,
        min_electron_temperature = cathode_electron_temp_eV,
        transition_function = LinearTransition(inner_outer_transition_length_m, 0.0),
        electron_pressure_coupled = false,
        scheme = HyperbolicScheme(;
            flux_function = eval(Symbol(flux_function)), limiter = eval(Symbol(limiter)), reconstruct
        ),
        solve_background_neutrals = solve_background_neutrals,
        background_pressure = background_pressure_Torr * u"Torr",
        background_neutral_temperature = background_temperature_K * u"K",
    )

    solution = run_simulation(
        config; ncells = num_cells, nsave = num_save,
        duration = duration_s, dt = dt_s, verbose = verbose, adaptive
    )

    return solution
end

function run_simulation(json_path::String; single_section = false, is_path = true, nonstandard_keys = false, verbose = true)

    if is_path
        json_content = JSON3.read(read(json_path, String))
    else
        json_content = JSON3.read(json_path)
    end

    run_simulation(json_content; single_section, nonstandard_keys, verbose)
end
