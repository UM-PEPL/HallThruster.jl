
abstract type InitialCondition end

struct DefaultInitialization <: InitialCondition end

#=
function IC!(U, z, fluids, L) #for testing light solve, energy equ is in eV*number*density
    mi = fluids[1].species.element.m
    un = 150.0
    #ρn = 5e-6/0.004/abs(un) - z / L * 5e-6/0.004/abs(un)
    ρn0 = 5e-6/0.004/abs(un)
    z1 = L/7
    z2 = L/2.5
    ρn = if z < z1
        ρn0
    elseif z < z2
        ρn0 - (z - z1) * (0.999 * ρn0) / (z2 - z1)
    else
        ρn0 / 1000
    end
    ui = if z < L/2 
        -1000 + 80000(z/L)^2
    else
        15000 + 10000z/L
    end

    ρi = mi * (2e17 + 9e17 * exp(-(4 * (z - L/4) / 0.033)^2))
    Tev = 3 + 37 * exp(-(2 * (z - L/2) / 0.023)^2)
    ne = ρi / mi

    ncharge = maximum(f.species.Z for f in fluids)

    if ncharge == 1
        U .= [ρn, ρi, ρi*ui, ne*Tev]
    elseif ncharge == 2
        U .= [ρn, ρi, ρi*ui, ρi/100, ρi/100 * sqrt(2) * ui, ne * Tev]
    elseif ncharge == 3
        U .= [ρn, ρi, ρi*ui, ρi/100, ρi/100 * sqrt(2) * ui, ρi/100, ρi/100 * sqrt(3) * ui, ne * Tev]
    end
    return U
end=#

initialize!(U, params) = initialize!(U, params, params.config.initial_condition)

function initialize!(U, params, model::InitialCondition)
    throw(ArgumentError("Function HallThruster.initialize!(U, params, model::$(typeof(model)) not yet implemented. For InitialCondition types other than DefaultInitialization(), this must be defined by the user!"))
end


function initialize!(U, params, ::DefaultInitialization)
    (;z_cell, config, index) = params
    (;ncharge, anode_Te, cathode_Te, domain, thruster, propellant, anode_mass_flow_rate, discharge_voltage) = config
    mi = propellant.m
    L_ch = thruster.geometry.channel_length

    ni_center = L_ch / 2
    ni_width = L_ch / 3
    ni_min = 2e17
    ni_max = 1.1e18
    ion_density_function = (z, Z) -> mi * (ni_min + (ni_max - ni_min) * exp(-((z - ni_center) / ni_width)^2)) / Z^2

    bohm_velocity = Z -> -sqrt(Z * e * 2/3 * anode_Te / mi)

    final_velocity = Z -> sqrt(2 * Z * e * discharge_voltage / mi)
    scale(Z) = 2/3 * (final_velocity(Z) - bohm_velocity(Z))
    ion_velocity_f1(z, Z) = bohm_velocity(Z) + scale(Z) * (z / L_ch)^2
    ion_velocity_f2(z, Z) = lerp(z, L_ch, domain[2], ion_velocity_f1(L_ch, Z), final_velocity(Z))

    ion_velocity_function = (z, Z) -> if z < L_ch
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

    neutral_function = z -> SmoothIf(transition_length = L_ch / 6)(z, L_ch / 2, ρn_0, ρn_1)

    number_density_function = z -> sum(Z * ion_density_function(z, Z) / mi for Z in 1:ncharge)

    Te_baseline = z -> lerp(z, domain[1], domain[2], anode_Te, cathode_Te)
    Te_min = min(anode_Te, cathode_Te)
    Te_max = config.discharge_voltage / 10
    Te_width = L_ch/3

    energy_function = z -> Te_baseline(z) + (Te_max - Te_min) * exp(-((z - L_ch) / Te_width)^2)

    # Fill the state vector
    for (i, z) in enumerate(z_cell)

        U[index.ρn, i] = neutral_function(z)

        for Z in 1:params.config.ncharge
            U[index.ρi[Z], i] = ion_density_function(z, Z)
            U[index.ρiui[Z], i] = ion_density_function(z, Z) * ion_velocity_function(z, Z)
        end

        U[index.nϵ, i] = number_density_function(z) * energy_function(z)
    end

    return U
end