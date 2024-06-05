function get_key(json_content, key, default)
    if haskey(json_content, key)
        return json_content[key]
    else
        return default
    end
end

function run_simulation(json_content::JSON3.Object; verbose = true)
    pressure_z0 = NaN
    pressure_dz = NaN
    pressure_pstar = NaN
    pressure_alpha = NaN
    apply_thrust_divergence_correction = true
    adaptive = true

    (;
        # Design
        channel_length, inner_radius, outer_radius,
        magnetic_field_file, magnetically_shielded,
        propellant, wall_material,
        anode_potential, cathode_potential,
        anode_mass_flow_rate,
        # Simulation
        anom_model, cathode_location_m, max_charge,
        num_cells, dt_s, duration_s, num_save,
        flux_function, limiter, reconstruct,
        ion_wall_losses, electron_ion_collisions,
        # Parameters
        sheath_loss_coefficient,
        ion_temp_K, neutral_temp_K, neutral_velocity_m_s,
        cathode_electron_temp_eV, inner_outer_transition_length_m,
        background_pressure_Torr, background_temperature_K,
    ) = json_content

    # Handle optional keys

    # Thruster name
    thruster_name = get_key(json_content, :thruster_name, "Unnamed Thruster")

    # Adaptive timestepping
    adaptive = get_key(json_content, :adaptive, true)

    # Anom coefficients
    anom_model_coeffs = get_key(json_content, :anom_model_coeffs, [0.00625, 0.0625])

    # Whether inbuilt divergence model is used to correct thrust
    apply_thrust_divergence_correction = get_key(json_content, :apply_thrust_divergence_correction, true)

    # neutral ingestion multiplier
    neutral_ingestion_multiplier = get_key(json_content, :neutral_ingestion_multiplier, 1.0)

    # Optional parameters for pressure-dependent models
    if  anom_model == "ShiftedTwoZone" || anom_model == "ShiftedTwoZoneBohm" ||
        anom_model == "ShiftedMultiBohm" || anom_model == "ShiftedGaussianBohm"
        (;pressure_z0, pressure_dz, pressure_pstar, pressure_alpha) = json_content
    end

    geometry = Geometry1D(;channel_length, outer_radius, inner_radius)

    if thruster_name == "SPT-100"
        bfield_func = HallThruster.B_field_SPT_100 $ (0.016, channel_length)
    else
        bfield_data = readdlm(magnetic_field_file, ',')
        bfield_func = HallThruster.LinearInterpolation(bfield_data[:, 1], bfield_data[:, 2])
    end

    thruster = HallThruster.Thruster(;
        name = thruster_name,
        geometry = geometry,
        magnetic_field = bfield_func,
        shielded = magnetically_shielded
    )

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
        coeff_tuple = (anom_model_coeffs[1], anom_model_coeffs[2])
        ShiftedTwoZoneBohm(coeff_tuple, pressure_z0, pressure_dz, pressure_pstar, pressure_alpha)
    elseif anom_model == "ShiftedMultiBohm"
        N = length(anom_model_coeffs)
        ShiftedMultiBohm(
            anom_model_coeffs[1:N÷2], anom_model_coeffs[N÷2 + 1], pressure_z0, pressure_dz, pressure_pstar, pressure_alpha
        )
    elseif anom_model == "ShiftedGaussianBohm"
        ShiftedGaussianBohm(
            anom_model_coeffs[1], anom_model_coeffs[2], anom_model_coeffs[3], anom_model_coeffs[4],
            pressure_z0, pressure_dz, pressure_pstar, pressure_alpha
        )
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
        transition_length = inner_outer_transition_length_m,
        scheme = HyperbolicScheme(;
            flux_function = eval(Symbol(flux_function)), limiter = eval(Symbol(limiter)), reconstruct
        ),
        background_pressure = background_pressure_Torr * u"Torr",
        background_neutral_temperature = background_temperature_K * u"K",
        neutral_ingestion_multiplier,
        apply_thrust_divergence_correction,
    )

    solution = run_simulation(
        config; ncells = num_cells, nsave = num_save,
        duration = duration_s, dt = dt_s, verbose = verbose, adaptive
    )

    return solution
end

function run_simulation(json_path::String; is_path = true, kwargs...)

    if is_path
        json_content = JSON3.read(read(json_path, String))
    else
        json_content = JSON3.read(json_path)
    end

    run_simulation(json_content; kwargs...)
end
