
abstract type InitialCondition end

struct DefaultInitialization <: InitialCondition end

initialize!(U, params) = initialize!(U, params, params.config.initial_condition)

function initialize!(U, params, model::InitialCondition)
    throw(ArgumentError("Function HallThruster.initialize!(U, params, model::$(typeof(model)) not yet implemented. For InitialCondition types other than DefaultInitialization(), this must be defined by the user!"))
end

function initialize!(U, params, ::DefaultInitialization)
    (;z_cell, config, index) = params
    (;ncharge, anode_Te, cathode_Te, domain, thruster, propellant, discharge_voltage, anode_mass_flow_rate) = config
    mi = propellant.m
    L_ch = thruster.geometry.channel_length
    z0 = domain[1]

    ni_center = L_ch / 2
    ni_width = L_ch / 3
    ni_min = 2e17
    ni_max = 1e18
    scaling_factor = sqrt(discharge_voltage / 300) * (anode_mass_flow_rate / 5e-6)
    ion_density_function = (z, Z) -> mi * scaling_factor * (ni_min + (ni_max - ni_min) * exp(-(((z-z0) - ni_center) / ni_width)^2)) / Z^2

    bohm_velocity = Z -> -sqrt(Z * e * anode_Te / mi)

    final_velocity = Z -> sqrt(2 * Z * e * discharge_voltage / mi)
    scale(Z) = 2/3 * (final_velocity(Z) - bohm_velocity(Z))
    ion_velocity_f1(z, Z) = bohm_velocity(Z) + scale(Z) * ((z-z0) / L_ch)^2
    ion_velocity_f2(z, Z) = lerp(z, z0 + L_ch, domain[2], ion_velocity_f1(L_ch, Z), final_velocity(Z))

    ion_velocity_function = (z, Z) -> if z - z0 < L_ch
        ion_velocity_f1(z, Z)
    else
        ion_velocity_f2(z, Z)
    end

    ρn_0 = inlet_neutral_density(config)
    # add recombined neutrals
    for Z in 1:config.ncharge
        ρn_0 -= ion_velocity_function(0.0, Z) * ion_density_function(0.0, Z) / config.neutral_velocity
    end

    ρn_1 = 0.01 * ρn_0

    neutral_function = z -> SmoothIf(transition_length = L_ch / 6)(z-z0, L_ch / 2, ρn_0, ρn_1)

    number_density_function = z -> sum(Z * ion_density_function(z, Z) / mi for Z in 1:ncharge)

    Te_baseline = z -> lerp(z, domain[1], domain[2], anode_Te, cathode_Te)
    Te_min = min(anode_Te, cathode_Te)
    Te_max = (config.discharge_voltage / 10)
    Te_width = L_ch/3

    energy_function = z -> 3/2 * (Te_baseline(z) + (Te_max - Te_min) * exp(-(((z-z0) - L_ch) / Te_width)^2))

    # Fill the state vector
    for (i, z) in enumerate(z_cell)

        U[index.ρn, i] = neutral_function(z)

        for Z in 1:params.config.ncharge
            U[index.ρi[Z], i] = ion_density_function(z, Z)
            U[index.ρiui[Z], i] = ion_density_function(z, Z) * ion_velocity_function(z, Z)
        end

        U[index.nϵ, i] = number_density_function(z) * energy_function(z)
    end

    # Initialize the anomalous collision frequency
    initialize_anom!(params.cache.νan, params.config.anom_model, U, params)

    return U
end