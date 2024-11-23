function get_key(json_content, key, default)
    if haskey(json_content, key)
        return json_content[key]
    else
        return default
    end
end

function config_from_json(json_content::JSON3.Object, dir = "")
    pressure_z0 = NaN
    pressure_dz = NaN
    pressure_pstar = NaN
    pressure_alpha = NaN
    apply_thrust_divergence_correction = true

    (;     # Design
    channel_length, inner_radius, outer_radius,
    magnetic_field_file, magnetically_shielded,
    propellant, wall_material,
    anode_potential, cathode_potential,
    anode_mass_flow_rate,     # Simulation
    anom_model, cathode_location_m, max_charge,
    flux_function, limiter, reconstruct,
    ion_wall_losses, electron_ion_collisions,     # Parameters
    sheath_loss_coefficient,
    ion_temp_K, neutral_temp_K, neutral_velocity_m_s,
    cathode_electron_temp_eV, inner_outer_transition_length_m,
    background_pressure_Torr, background_temperature_K) = json_content

    # Thruster name
    thruster_name = get_key(json_content, :thruster_name, "Unnamed Thruster")

    # Anom coefficients
    anom_model_coeffs = get_key(json_content, :anom_model_coeffs, [0.00625, 0.0625])

    # Whether inbuilt divergence model is used to correct thrust
    apply_thrust_divergence_correction = get_key(json_content,
        :apply_thrust_divergence_correction, true,)

    # Whether we solve a quasi-1D plume expansion
    solve_plume = get_key(json_content, :solve_plume, true)

    # Electron loss coeffient in the plume
    electron_plume_loss_scale = get_key(json_content, :plume_loss_coefficient, 1.0)

    # neutral ingestion multiplier
    neutral_ingestion_multiplier::Float64 = get_key(json_content,
        :neutral_ingestion_multiplier, 1.0,)

    # Optional parameters for pressure-dependent models
    if anom_model == "ShiftedTwoZone" || anom_model == "ShiftedTwoZoneBohm" ||
       anom_model == "ShiftedMultiBohm" || anom_model == "ShiftedGaussianBohm"
        (; pressure_z0, pressure_dz, pressure_pstar, pressure_alpha) = json_content
    end

    geometry = Geometry1D(; channel_length, outer_radius, inner_radius)

    # Look for magnetic field file in same dir as JSON file if not found in dir where script was run
    include_paths = ["", dir]
    for path in include_paths
        bfield_file = joinpath(path, magnetic_field_file)
        if ispath(bfield_file)
            magnetic_field_file = bfield_file
            break
        end
    end

    bfield_func = try
        MagneticField(magnetic_field_file)
    catch e
        if thruster_name == "SPT-100"
            @warn "Could not find magnetic field file at path $(magnetic_field_file). \
                 Using default SPT-100 field."
            SPT_100.magnetic_field
        else
            error(e)
        end
    end

    thruster = HallThruster.Thruster(;
        name = thruster_name,
        geometry = geometry,
        magnetic_field = bfield_func,
        shielded = magnetically_shielded,)

    # assign anomalous transport model
    anom_model = if anom_model == "NoAnom"
        NoAnom()
    elseif anom_model == "Bohm"
        Bohm(anom_model_coeffs[1])
    elseif anom_model == "TwoZoneBohm"
        TwoZoneBohm((anom_model_coeffs[1], anom_model_coeffs[2]))
    elseif anom_model == "MultiLogBohm"
        MultiLogBohm(anom_model_coeffs)
    elseif anom_model == "ShiftedTwoZone" || anom_model == "ShiftedTwoZoneBohm"
        c1, c2 = anom_model_coeffs[1], anom_model_coeffs[2]
        LogisticPressureShift(
            TwoZoneBohm(c1, c2), pressure_z0, pressure_dz, pressure_pstar, float(pressure_alpha),)
    elseif anom_model == "ShiftedMultiBohm"
        N = length(anom_model_coeffs)
        LogisticPressureShift(
            MultiLogBohm(anom_model_coeffs[1:(N ÷ 2)], anom_model_coeffs[N ÷ 2 + 1]),
            pressure_z0, pressure_dz, pressure_pstar, float(pressure_alpha),)
    elseif anom_model == "ShiftedGaussianBohm"
        LogisticPressureShift(
            GaussianBohm(anom_model_coeffs[1], anom_model_coeffs[2],
                anom_model_coeffs[3], anom_model_coeffs[4],),
            pressure_z0, pressure_dz, pressure_pstar, float(pressure_alpha),)
    end

    config = HallThruster.Config(;
        thruster,
        propellant = eval(Symbol(propellant)),
        anom_model,
        domain = (0.0, cathode_location_m),
        discharge_voltage = anode_potential,
        anode_mass_flow_rate = anode_mass_flow_rate,
        cathode_potential = cathode_potential,
        ncharge = max_charge,
        wall_loss_model = WallSheath(eval(Symbol(wall_material)),
            Float64(sheath_loss_coefficient),),
        ion_wall_losses = ion_wall_losses,
        anode_Tev = cathode_electron_temp_eV,
        cathode_Tev = cathode_electron_temp_eV,
        LANDMARK = false,
        ion_temperature_K = ion_temp_K,
        neutral_temperature_K = neutral_temp_K,
        neutral_velocity = neutral_velocity_m_s,
        electron_ion_collisions = electron_ion_collisions,
        transition_length = inner_outer_transition_length_m,
        scheme = HyperbolicScheme(;
            flux_function = eval(Symbol(flux_function)),
            limiter = eval(Symbol(limiter)),
            reconstruct,),
        background_pressure_Torr = background_pressure_Torr,
        background_temperature_K = background_temperature_K,
        neutral_ingestion_multiplier,
        solve_plume,
        apply_thrust_divergence_correction,
        electron_plume_loss_scale,)

    return config
end

function run_simulation(json_content::JSON3.Object, dir::String = ""; verbose = true)
    adaptive = get_key(json_content, :adaptive, true)
    num_save = get_key(json_content, :num_save, 100)
    num_cells = get_key(json_content, :num_cells, 200)
    duration_s = get_key(json_content, :duration_s, 1e-3)
    dt_s = get_key(json_content, :dt_s, 1e-8)
    dtmin = get_key(json_content, :min_dt_s, 1e-10)
    dtmax = get_key(json_content, :max_dt_s, 1e-7)
    max_small_steps = get_key(json_content, :max_small_steps, 100)

    config = config_from_json(json_content, dir)

    solution = run_simulation(config; grid = EvenGrid(num_cells), nsave = num_save,
        duration = duration_s, dt = dt_s, verbose = verbose, adaptive,
        dtmin, dtmax, max_small_steps,)

    return solution
end

function run_simulation(json_path::String; is_path = true, kwargs...)
    if is_path
        json_content = JSON3.read(read(json_path, String))
    else
        json_content = JSON3.read(json_path)
    end

    return run_simulation(json_content, dirname(abspath(json_path)); kwargs...)
end
