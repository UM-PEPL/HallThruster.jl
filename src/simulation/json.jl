function frame_dict(sol, frame)
    filler = zeros(length(sol.params.z_cell))
    ncharge = sol.params.config.ncharge
    return Dict(
        "thrust" => thrust(sol, frame),
        "discharge_current" => discharge_current(sol, frame),
        "mass_eff" => mass_eff(sol, frame),
        "voltage_eff" => voltage_eff(sol, frame),
        "current_eff" => current_eff(sol, frame),
        "divergence_eff" => divergence_eff(sol, frame),
        "t" => sol.t[frame],
        "z" => sol.params.z_cell,
        "nn" => sol[:nn, 1][frame],
        "ni_1" => sol[:ni, 1][frame],
        "ni_2" => (ncharge > 1 ? sol[:ni, 2][frame] : filler),
        "ni_3" => (ncharge > 2 ? sol[:ni, 3][frame] : filler),
        "ne" => sol[:ne][frame],
        "ui_1" => sol[:ui, 1][frame],
        "ui_2" => (ncharge > 1 ? sol[:ui, 2][frame] : filler),
        "ui_3" => (ncharge > 2 ? sol[:ui, 3][frame] : filler),
        "niui_1" => sol[:niui, 1][frame],
        "niui_2" => (ncharge > 1 ? sol[:niui, 2][frame] : filler),
        "niui_3" => (ncharge > 2 ? sol[:niui, 3][frame] : filler),
        "ue" => sol[:ue][frame],
        "potential" => sol[:potential][frame],
        "electric_field" => sol[:electric_field][frame],
        "Tev" => sol[:Tev][frame],
        "pe" => sol[:pe][frame],
        "grad_pe" => sol[:grad_pe][frame],
        "nu_en" => sol[:νen][frame],
        "nu_ei" => sol[:νei][frame],
        "nu_anom" => sol[:νan][frame],
        "nu_class" => sol[:νc][frame],
        "mobility" => sol[:mobility][frame]
    )
end

function write_to_json(filename, sol)
    num_frames = length(sol.t)
    frames = map(frame -> frame_dict(sol, frame), 1:num_frames)
    return JSON3.write(filename, frames)
end

function get_key(json_content, key, default)
    if haskey(json_content, key)
        return json_content[key]
    else
        return default
    end
end

function assign_json_key(keyname, opt, pairs...)
    sym = Symbol(opt)
    dict = Dict(pairs...)
    if !haskey(dict, sym)
        valid_options = map(string, [k for (k, _) in dict])
        throw(ArgumentError("Invalid $(keyname) '$(opt)' for JSON input. \
                            Please select one of $(valid_options)"))
    end
    return dict[sym]
end

function config_from_json(json_content::JSON3.Object; verbose = true)
    (;     # Design
    channel_length,
    inner_radius,
    outer_radius,
    magnetic_field_file,
    magnetically_shielded,
    propellant,
    wall_material,
    anode_potential,
    cathode_potential,
    anode_mass_flow_rate,     # Simulation
    anom_model,
    cathode_location_m,
    max_charge,
    flux_function,
    limiter,
    reconstruct,
    ion_wall_losses,
    electron_ion_collisions,     # Parameters
    sheath_loss_coefficient,
    ion_temp_K,
    neutral_temp_K,
    neutral_velocity_m_s,
    cathode_electron_temp_eV,
    inner_outer_transition_length_m,
    background_pressure_Torr,
    background_temperature_K
) = json_content

    # Thruster name
    thruster_name = get_key(json_content, :thruster_name, "Unnamed Thruster")

    # Anom coefficients
    anom_model_coeffs = get_key(json_content, :anom_model_coeffs, [0.00625, 0.0625])

    # Whether inbuilt divergence model is used to correct thrust
    apply_thrust_divergence_correction = get_key(
        json_content, :apply_thrust_divergence_correction, true
    )

    # Whether we solve a quasi-1D plume expansion
    solve_plume = get_key(json_content, :solve_plume, true)

    # Electron loss coeffient in the plume
    electron_plume_loss_scale = get_key(json_content, :plume_loss_coefficient, 1.0)

    # neutral ingestion multiplier
    neutral_ingestion_multiplier::Float64 = get_key(
        json_content, :neutral_ingestion_multiplier, 1.0
    )

    geometry = Geometry1D(; channel_length, outer_radius, inner_radius)

    bfield_func = try
        MagneticField(magnetic_field_file)
    catch e
        if thruster_name == "SPT-100"
            if verbose
                @warn "Could not find provided magnetic field file. Using default SPT-100 field."
            end
            B_field_SPT_100()
        else
            error(e)
        end
    end

    # Construct a `thruster`
    thruster = HallThruster.Thruster(;
        name = thruster_name,
        geometry = geometry,
        magnetic_field = bfield_func,
        shielded = magnetically_shielded
    )

    # Assign anomalous transport models
    anom_model_fn = assign_json_key(
        "anom_model",
        anom_model,
        :NoAnom => (_, _) -> NoAnom(),
        :Bohm => (_, c) -> Bohm(c[1]),
        :TwoZoneBohm => (_, c) -> TwoZoneBohm(c[1], c[2]),
        :MultiLogBohm => (_, c) -> MultiLogBohm(c),
        :ShiftedTwoZone => (js, c) -> ShiftedTwoZoneBohm(
            (c[1], c[2]),
            js.pressure_z0,
            js.pressure_dz,
            js.pressure_pstar,
            js.pressure_alpha
        ),
        :ShiftedMultiBohm => (js, c) -> ShiftedMultiBohm(
            c[1:(length(c) ÷ 2)],
            c[length(c) ÷ 2 + 1],
            js.pressure_z0,
            js.pressure_dz,
            js.pressure_pstar,
            js.pressure_alpha
        ),
        :ShiftedGaussianBohm => (js, c) -> ShiftedGaussianBohm(
            c[1],
            c[2],
            c[3],
            c[4],
            js.pressure_z0,
            js.pressure_dz,
            js.pressure_pstar,
            js.pressure_alpha
        )
    )
    anom_model = anom_model_fn(json_content, anom_model_coeffs)

    propellant = assign_json_key(
        "propellant", propellant, :Xenon => Xenon, :Krypton => Krypton, :Argon => Argon
    )

    wall_material = assign_json_key(
        "wall material", wall_material, :BNSiO2 => BNSiO2, :BoronNitride => BoronNitride
    )

    flux_function = assign_json_key(
        "flux function",
        flux_function,
        :rusanov => rusanov,
        :HLLE => HLLE,
        :global_lax_friedrichs => global_lax_friedrichs
    )

    limiter = assign_json_key("limiter", limiter, :van_leer => van_leer, :minmod => minmod)

    config = HallThruster.Config(;
        thruster,
        propellant = propellant,
        anom_model,
        domain = (0.0, cathode_location_m),
        discharge_voltage = anode_potential - cathode_potential,
        anode_mass_flow_rate = anode_mass_flow_rate,
        cathode_potential = cathode_potential,
        ncharge = max_charge,
        wall_loss_model = WallSheath(wall_material, Float64(sheath_loss_coefficient)),
        ion_wall_losses = ion_wall_losses,
        cathode_Te = cathode_electron_temp_eV,
        LANDMARK = false,
        ion_temperature = ion_temp_K,
        neutral_temperature = neutral_temp_K,
        neutral_velocity = neutral_velocity_m_s,
        electron_ion_collisions = electron_ion_collisions,
        min_electron_temperature = cathode_electron_temp_eV,
        transition_length = inner_outer_transition_length_m,
        scheme = HyperbolicScheme(; flux_function, limiter, reconstruct),
        background_pressure = background_pressure_Torr,
        background_neutral_temperature = background_temperature_K,
        neutral_ingestion_multiplier,
        solve_plume,
        apply_thrust_divergence_correction,
        electron_plume_loss_scale
    )

    return config
end

function config_from_json(json_path::String; is_path = true, kwargs...)
    if is_path
        json_content = JSON3.read(read(json_path, String))
    else
        json_content = JSON3.read(json_path)
    end

    return config_from_json(json_content)
end

function run_simulation(json_content::JSON3.Object; verbose = true)
    adaptive = get_key(json_content, :adaptive, true)
    num_save = get_key(json_content, :num_save, 100)
    num_cells = get_key(json_content, :num_cells, 200)
    duration_s = get_key(json_content, :duration_s, 1e-3)
    dt_s = get_key(json_content, :dt_s, 1e-8)
    dtmin = get_key(json_content, :min_dt_s, 1e-10)
    dtmax = get_key(json_content, :max_dt_s, 1e-7)
    max_small_steps = get_key(json_content, :max_small_steps, 100)

    config = config_from_json(json_content; verbose)

    solution = run_simulation(
        config;
        grid = EvenGrid(num_cells),
        nsave = num_save,
        duration = duration_s,
        dt = dt_s,
        verbose = verbose,
        adaptive,
        dtmin,
        dtmax,
        max_small_steps
    )

    return solution
end

function run_simulation(json_path::String; is_path = true, kwargs...)
    if is_path
        json_content = JSON3.read(read(json_path, String))
    else
        json_content = JSON3.read(json_path)
    end

    return run_simulation(json_content; kwargs...)
end
