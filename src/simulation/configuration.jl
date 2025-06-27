"""
$(TYPEDEF)
Hall thruster configuration struct. Only four mandatory fields: `discharge_voltage`, `thruster`, `anode_mass_flow_rate`, and `domain`.

## Mandatory Fields
$(TYPEDFIELDS)
"""
struct Config{A <: AnomalousTransportModel, TC <: ThermalConductivityModel, W <: WallLossModel, IC <: InitialCondition, S_HS, S_E}
    """
    The thruster to simulate. See [Thrusters](@ref) for more information
    """
    thruster::Thruster
    """
    The simulation domain, given as (left, right). Dimensions are in meters.
    """
    domain::Tuple{Float64, Float64}
    """
    The potential difference between anode and cathode, in V. Used to set the left boundary condition for electrostatic potential.
    """
    discharge_voltage::Float64
    """
    The propellants to be used. See [Propellants](propellants.md) for more.
    """
    propellants::Vector{Propellant}
    """
    ---
    # Optional fields
    ---
    The potential at the right boundary of the simulation. **Default:** `0`
    """
    cathode_coupling_voltage::Float64
    """
    Can be either `:sheath` or `:dirichlet`. If `:sheath`, electron temperature has a Neumann boundary condition at the anode and a self-consistent anode sheath potential is computed. If `:dirichlet`, electron temperature at anode is set to `anode_Tev` and no sheath potential is modeled. **Default:** `:sheath`.
    """
    anode_boundary_condition::Symbol
    """
    Electron temperature at left boundary (anode) if `anode_boundary_condition == :sheath`. **Default:** 2.0.
    """
    anode_Tev::Float64
    """
    Electron temperature at right boundary (cathode). **Default:** 2.0
    """
    cathode_Tev::Float64
    """
    Model for computing the anomalous collision frequency. See [Anomalous Transport](../reference/anomalous_transport.md) for more info. **Default:** `TwoZoneBohm(1/160, 1/16)`.
    """
    anom_model::A
    """
    How radial losses due to sheaths are computed. Other wall models are described on the [Wall Loss Models](wall_loss_models.md) page. **Default:** `WallSheath(BNSiO2, 1.0)`.
    """
    wall_loss_model::W
    """
    Model for the cross-field electron thermal conductivity. See [Electron Thermal Conductivity](../reference/electron_thermal_conductivity.md) for more. **Default:** `Mitchner()`
    """
    conductivity_model::TC
    """
    Whether to include electron-ion collisions. See [Collisions and Reactions](@ref) for more. **Default:** `true`.
    """
    electron_ion_collisions::Bool
    """
    Whether we model ion losses to the walls. **Default:** `false`.
    """
    ion_wall_losses::Bool
    """
    The pressure of the background neutrals, in Pascals. These background neutrals are injected at the anode to simulate the ingestion of facility neutrals. **Default:** `0.0`
    """
    background_pressure_Torr::Float64
    """
    The temperature of the background neutrals, in K. **Default:** `150.0`.
    """
    background_temperature_K::Float64
    """
    The factor by which the ingested mass flow rate computed from the background pressure and temperature is multiplied. **Default:** 1.
    """
    neutral_ingestion_multiplier::Float64
    """
    Whether quasi-1D beam expansion should be modelled outside of the channel. See [Quasi-1D plume model](../explanation/plume.md) for more. **Default:** `false`
    """
    solve_plume::Bool
    """
    Whether the thrust output by HallThruster.jl should include a divergence correction factor of `cos(divergence_angle)`. **Default:** `false`.
    """
    apply_thrust_divergence_correction::Bool
    """
    The degree to which radial electron losses are applied in the plume. See [Wall Loss Models](@ref) for more information. **Default:** 1.
    """
    electron_plume_loss_scale::Float64
    """
    Factor by which the magnetic field is increased or decreased compared to the one in the provided `Thruster` struct. **Default:** `1.0`.
    """
    magnetic_field_scale::Float64
    """
    Distance over which the transition between inside and outside the channel is smoothed. Affects wall losses as well as two-zone Bohm-like transport models. **Default:** `0.1 * thruster.geometry.channel_length`
    """
    transition_length::Float64
    """
    Whether to employ gradient reconstruction
    """
    reconstruct::Bool
    """
    An `InitialCondition`; see [Initialization](../explanation/initialization.md) for more information. **Default:** `DefaultInitialization()`.
    """
    initial_condition::IC
    """
    The degree to which the energy is solved implicitly. `0.0` is a fully-explicit forward Euler, `0.5` is Crank-Nicholson, and `1.0` is backward Euler. **Default:** `1.0`.
    """
    implicit_energy::Float64
    """
    Additional directories in which we should look for rate coefficients. These are searched before the default directory, so replacement rate coefficients for built-in propellants can be provided. **Default:** String[]
    """
    reaction_rate_directories::Vector{String}
    """
     How many times to smooth the anomalous transport profile. Only useful for transport models that depend on the plasma properties. **Default:** `0`

    ---
    # Verification and validation options
    ---
    These options are used in code benchmarking and verification and are not usually needed by end-users.
    See [Verification and validation](../explanation/verification.md) for an explanation of how we use these to verify the accuracy of the code.
    """
    anom_smoothing_iters::Int
    """
    Whether we are using the physics model from the LANDMARK benchmark. This affects whether certain terms are included in the equations, such as electron and heavy species momentum transfer due to ionization and the form of the electron thermal conductivity. **Default:** `false`.
    """
    LANDMARK::Bool
    """
    Model for ionization reactions. **Default:** `:Lookup`.
    """
    ionization_model::Symbol
    """
    Model for excitation reactions. **Default:** `:Lookup`. 
    """
    excitation_model::Symbol
    """
    Model for elastic scattering collisions between electrons and neutral atoms. **Default:** `:Lookup`.
    """
    electron_neutral_model::Symbol
    """
    Extra source term for heavy species. **Default:** `nothing`
    """
    source_heavy_species::S_HS
    """
    Extra source term for electron energy equation. **Default:** `nothing`. 
    """
    source_energy::S_E

    function Config(;
            # Mandatory arguments
            thruster::Thruster,
            domain,
            discharge_voltage,
            propellants = nothing,
            # Optional arguments
            cathode_coupling_voltage = 0.0,
            anode_boundary_condition = :sheath,
            cathode_Tev = 2.0,
            anode_Tev = cathode_Tev,
            anom_model::A = TwoZoneBohm(1 / 160, 1 / 16),
            wall_loss_model::W = WallSheath(BNSiO2, 1.0),
            conductivity_model::TC = Mitchner(),
            electron_ion_collisions = true,
            ion_wall_losses = false,
            background_pressure_Torr = 0.0,
            background_temperature_K = 100.0,
            neutral_ingestion_multiplier = 1.0,
            solve_plume = false,
            apply_thrust_divergence_correction = false,
            electron_plume_loss_scale = 1.0,
            magnetic_field_scale::Float64 = 1.0,
            transition_length = 0.1 * thruster.geometry.channel_length,
            reconstruct::Bool = true,
            initial_condition::IC = DefaultInitialization(),
            implicit_energy = 1.0,
            reaction_rate_directories = String[],
            anom_smoothing_iters = 0,
            LANDMARK = false,
            ionization_model = :Lookup,
            excitation_model = :Lookup,
            electron_neutral_model = :Lookup,
            source_heavy_species = Returns(0.0),
            source_energy = Returns(0.0),
            # Backwards-compatible arguments
            anode_mass_flow_rate = nothing,
            propellant = Xenon,
            neutral_temperature_K = nothing,
            neutral_velocity = nothing,
            ion_temperature_K = DEFAULT_ION_TEMPERATURE_K,
            ncharge = 1,
        ) where {
            A <: AnomalousTransportModel,
            TC <: ThermalConductivityModel,
            W <: WallLossModel,
            IC <: InitialCondition,
        }

        # Set up propellants
        if isnothing(propellants)
            if isnothing(anode_mass_flow_rate)
                error("Must supply either a vector of propellants or an anode mass flow rate")
            end
            prop = Propellant(
                propellant, anode_mass_flow_rate;
                max_charge = ncharge, velocity_m_s = neutral_velocity,
                temperature_K = neutral_temperature_K, ion_temperature_K
            )
            propellants = [prop]
        end

        # Convert to Float64 if using Unitful
        discharge_voltage = convert_to_float64(discharge_voltage, units(:V))
        cathode_coupling_voltage = convert_to_float64(cathode_coupling_voltage, units(:V))

        anode_Tev = convert_to_float64(anode_Tev, units(:eV))

        domain = (
            convert_to_float64(domain[1], units(:m)),
            convert_to_float64(domain[2], units(:m)),
        )

        background_temperature_K = convert_to_float64(background_temperature_K, units(:K))
        background_pressure_Torr = convert_to_float64(background_pressure_Torr, units(:Pa))

        transition_length = convert_to_float64(transition_length, units(:m))

        if anode_boundary_condition ∉ [:sheath, :dirichlet]
            throw(ArgumentError("Anode boundary condition must be one of [:sheath, :dirichlet]. Got: $(anode_boundary_condition)"))
        end

        return new{A, TC, W, IC, typeof(source_heavy_species), typeof(source_energy)}(
            # Mandatory arguments
            thruster,
            domain,
            discharge_voltage,
            propellants,
            # Optional arguments
            cathode_coupling_voltage,
            anode_boundary_condition,
            anode_Tev,
            cathode_Tev,
            anom_model,
            wall_loss_model,
            conductivity_model,
            electron_ion_collisions,
            ion_wall_losses,
            background_pressure_Torr,
            background_temperature_K,
            neutral_ingestion_multiplier,
            solve_plume,
            apply_thrust_divergence_correction,
            electron_plume_loss_scale,
            magnetic_field_scale,
            transition_length,
            reconstruct,
            initial_condition,
            implicit_energy,
            reaction_rate_directories,
            anom_smoothing_iters,
            LANDMARK,
            ionization_model,
            excitation_model,
            electron_neutral_model,
            source_heavy_species,
            source_energy,
        )
    end
end

function params_from_config(config)
    # TODO: make work better with mutliple propellants
    un_B = background_neutral_velocity(config)
    nn_B = background_neutral_density(config)
    un = config.propellants[1].velocity_m_s
    ingestion_density = nn_B * un_B / un * config.neutral_ingestion_multiplier

    return (;
        # Copied directly from config
        propellants = config.propellants,
        reconstruct = config.reconstruct,
        thruster = config.thruster,
        anode_bc = config.anode_boundary_condition,
        landmark = config.LANDMARK,
        transition_length = config.transition_length,
        Te_L = config.anode_Tev,
        Te_R = config.cathode_Tev,
        implicit_energy = config.implicit_energy,
        ingestion_density,
        ion_temperature_K = config.propellants[1].ion_temperature_K,
        neutral_velocity = un,
        anode_mass_flow_rate = config.propellants[1].flow_rate_kg_s,
        neutral_ingestion_multiplier = config.neutral_ingestion_multiplier,
        ion_wall_losses = config.ion_wall_losses,
        wall_loss_scale = wall_loss_scale(config.wall_loss_model),
        plume_loss_scale = config.electron_plume_loss_scale,
        anom_smoothing_iters = config.anom_smoothing_iters,
        discharge_voltage = config.discharge_voltage,
        cathode_coupling_voltage = config.cathode_coupling_voltage,
        electron_ion_collisions = config.electron_ion_collisions,
    )
end

#=============================================================================
 Serialization of Config to JSON
==============================================================================#

# Don't write source terms to output or read them from input
function Serialization.exclude(::Type{C}) where {C <: Config}
    return (:source_heavy_species, :source_energy)
end

function Serialization.deserialize(::Type{C}, x) where {C <: Config}
    d = copy(x)
    # Handle configs from older versions
    if !haskey(d, :propellants)
        gas = get(d, :propellant, Xenon)
        max_charge = get(d, :ncharge, 1)
        velocity_m_s = get(d, :neutral_velocity, nothing)
        temperature_K = get(d, :neutral_temperature_K, nothing)
        ion_temperature_K = get(d, :ion_temperature_K, DEFAULT_ION_TEMPERATURE_K)
        flow_rate_kg_s = d[:anode_mass_flow_rate]
        prop = Propellant(gas, flow_rate_kg_s; max_charge, velocity_m_s, temperature_K, ion_temperature_K)
        prop_dict = serialize(prop)
        for key in [
                :propellant, :ncharge, :anode_mass_flow_rate,
                :neutral_velocity, :neutral_temperature_K, :ion_temperature_K,
            ]
            delete!(d, key)
        end
        d[:propellants] = [prop_dict]
    end
    return deserialize(Serialization.Struct(), C, d)
end
