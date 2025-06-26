"""
$(TYPEDEF)
Hall thruster configuration struct. Only four mandatory fields: `discharge_voltage`, `thruster`, `anode_mass_flow_rate`, and `domain`.

## Mandatory Fields
$(TYPEDFIELDS)
"""
struct Config{A <: AnomalousTransportModel, TC <: ThermalConductivityModel, W <: WallLossModel, IC <: InitialCondition, S_N, S_IC, S_IM, S_E}
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
    The mass flow rate of neutral atoms through the anode, in kg/s.

    ---
    # Optional fields
    ---
    """
    anode_mass_flow_rate::Float64
    """
    Maximum ion charge state. **Default:** 1.
    """
    ncharge::Int
    """
    A `Gas`. See [Propellants](propellants.md) for more. **Default:** `Xenon`.
    """
    propellant::Gas
    """
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
    Neutral velocity in m/s. **Default:** `300.0`, or if `neutral_temperature` is set, that parameter is used to compute the velocity using a one-sided maxwellian flux approximation.
    """
    neutral_velocity::Float64
    """
    Neutral temperature in Kelvins. **Default:** `500.0`.
    """
    neutral_temperature_K::Float64
    """
    Ion temperature in Kelvins. **Default:** 1000.0
    """
    ion_temperature_K::Float64
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
    Extra user-provided neutral source term. **Default:** `nothing` 
    """
    source_neutrals::S_N
    """
    Vector of extra source terms for ion continuity, one for each charge state. **Default:** `nothing`.
    """
    source_ion_continuity::S_IC
    """
    Vector of extra source terms for ion momentum, one for each charge state. **Default:** `nothing`.
    """
    source_ion_momentum::S_IM
    """
    Extra source term for electron energy equation. **Default:** `nothing`. 
    """
    source_energy::S_E

    function Config(;
            # Mandatory arguments
            thruster::Thruster,
            domain,
            discharge_voltage,
            anode_mass_flow_rate,
            # Optional arguments
            ncharge = 1,
            propellant = Xenon,
            cathode_coupling_voltage = 0.0,
            anode_boundary_condition = :sheath,
            cathode_Tev = 2.0,
            anode_Tev = cathode_Tev,
            anom_model::A = TwoZoneBohm(1 / 160, 1 / 16),
            wall_loss_model::W = WallSheath(BNSiO2, 1.0),
            conductivity_model::TC = Mitchner(),
            electron_ion_collisions = true,
            neutral_velocity = nothing,
            neutral_temperature_K = nothing,
            ion_temperature_K = 1000.0,
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
            source_neutrals = nothing,
            source_ion_continuity = nothing,
            source_ion_momentum = nothing,
            source_energy = Returns(0.0),
        ) where {
            A <: AnomalousTransportModel,
            TC <: ThermalConductivityModel,
            W <: WallLossModel,
            IC <: InitialCondition,
        }
        # check that number of ion source terms matches number of charges for both
        # continuity and momentum
        #
        source_ion_continuity = ion_source_terms(
            ncharge, source_ion_continuity, "continuity",
        )
        source_ion_momentum = ion_source_terms(ncharge, source_ion_momentum, "momentum")

        # Neutral source terms
        if isnothing(source_neutrals)
            source_neutrals = Returns(0.0)
        end

        # Convert to Float64 if using Unitful
        discharge_voltage = convert_to_float64(discharge_voltage, units(:V))
        cathode_coupling_voltage = convert_to_float64(cathode_coupling_voltage, units(:V))

        anode_Tev = convert_to_float64(anode_Tev, units(:eV))
        cathode_Tev = convert_to_float64(cathode_Tev, units(:eV))

        default_neutral_velocity = 150.0 # m/s
        default_neutral_temp = 500.0 # K
        if isnothing(neutral_velocity) && isnothing(neutral_temperature_K)
            neutral_velocity = default_neutral_velocity
            neutral_temperature_K = default_neutral_temp
        elseif isnothing(neutral_temperature_K)
            neutral_temperature_K = default_neutral_temp
            neutral_velocity = convert_to_float64(neutral_velocity, units(:m) / units(:s))
        elseif isnothing(neutral_velocity)
            # compute neutral velocity from thermal speed
            neutral_temperature_K = convert_to_float64(neutral_temperature_K, units(:K))
            neutral_velocity = 0.25 *
                sqrt(8 * kB * neutral_temperature_K / π / propellant.m)
        else
            neutral_velocity = convert_to_float64(neutral_velocity, units(:m) / units(:s))
            neutral_temperature_K = convert_to_float64(neutral_temperature_K, units(:K))
        end

        ion_temperature_K = convert_to_float64(ion_temperature_K, units(:K))
        domain = (
            convert_to_float64(domain[1], units(:m)),
            convert_to_float64(domain[2], units(:m)),
        )

        anode_mass_flow_rate = convert_to_float64(
            anode_mass_flow_rate, units(:kg) / units(:s),
        )

        background_temperature_K = convert_to_float64(background_temperature_K, units(:K))
        background_pressure_Torr = convert_to_float64(background_pressure_Torr, units(:Pa))

        transition_length = convert_to_float64(transition_length, units(:m))

        if anode_boundary_condition ∉ [:sheath, :dirichlet]
            throw(ArgumentError("Anode boundary condition must be one of [:sheath, :dirichlet]. Got: $(anode_boundary_condition)"))
        end

        return new{
            A, TC, W, IC,
            typeof(source_neutrals), typeof(source_ion_continuity),
            typeof(source_ion_momentum), typeof(source_energy),
        }(
            # Mandatory arguments
            thruster,
            domain,
            discharge_voltage,
            anode_mass_flow_rate,
            # Optional arguments
            ncharge,
            propellant,
            cathode_coupling_voltage,
            anode_boundary_condition,
            anode_Tev,
            cathode_Tev,
            anom_model,
            wall_loss_model,
            conductivity_model,
            electron_ion_collisions,
            neutral_velocity,
            neutral_temperature_K,
            ion_temperature_K,
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
            source_neutrals,
            source_ion_continuity,
            source_ion_momentum,
            source_energy,
        )
    end
end

#=============================================================================
 Serialization of Config to JSON
==============================================================================#

# Don't write source terms to output or read them from input
function Serialization.exclude(::Type{C}) where {C <: Config}
    return (
        :source_neutrals, :source_ion_continuity,
        :source_ion_momentum, :source_potential, :source_energy,
    )
end

function ion_source_terms(ncharge, source, type)
    if ncharge != length(source)
        throw(ArgumentError("Number of ion $type source terms must match number of charges"))
    end
    return source
end

ion_source_terms(ncharge, ::Nothing, args...) = fill(Returns(0.0), ncharge)

function make_keys(fluid_range, subscript)
    len = length(fluid_range)
    if len == 1
        return (Symbol("ρ$(subscript)"))
    elseif len == 2
        return (
            Symbol("ρ$(subscript)"),
            Symbol("ρ$(subscript)u$(subscript)"),
        )
    elseif len == 3
        return (
            Symbol("ρ$(subscript)"),
            Symbol("ρ$(subscript)u$(subscript)"),
            Symbol("ρ$(subscript)E$(subscript)"),
        )
    else
        throw(ArgumentError("Too many equations on fluid (this should be unreachable)"))
    end
end

function configure_fluids(config)
    propellant = config.propellant

    neutral_fluid = ContinuityOnly(
        propellant(0); u = config.neutral_velocity, T = config.neutral_temperature_K,
    )
    ion_fluids = [
        IsothermalEuler(propellant(Z); T = config.ion_temperature_K)
            for Z in 1:(config.ncharge)
    ]

    fluids = [neutral_fluid; ion_fluids]

    species = [f.species for f in fluids]

    fluid_ranges = ranges(fluids)
    species_range_dict = Dict(Symbol(fluid.species) => 0:0 for fluid in fluids)

    for (fluid, fluid_range) in zip(fluids, fluid_ranges)
        species_range_dict[Symbol(fluid.species)] = fluid_range
    end

    last_fluid_index = fluid_ranges[end][end]
    is_velocity_index = fill(false, last_fluid_index)
    for i in 3:2:last_fluid_index
        is_velocity_index[i] = true
    end

    return fluids, fluid_ranges, species, species_range_dict, is_velocity_index
end

function configure_index(fluids, fluid_ranges)
    first_ion_fluid_index = findfirst(x -> x.species.Z > 0, fluids)

    keys_neutrals = (:ρn,)
    values_neutrals = (1,)

    keys_ions = (:ρi, :ρiui)
    values_ions = (
        [f[1] for f in fluid_ranges[first_ion_fluid_index:end]],
        [f[2] for f in fluid_ranges[first_ion_fluid_index:end]],
    )

    keys_fluids = (keys_neutrals..., keys_ions...)
    values_fluids = (values_neutrals..., values_ions...)
    index = NamedTuple{keys_fluids}(values_fluids)
    return index
end

function params_from_config(config)
    return (;
        # Copied directly from config
        thruster = config.thruster,
        ncharge = config.ncharge,
        mi = config.propellant.m,
        anode_bc = config.anode_boundary_condition,
        landmark = config.LANDMARK,
        transition_length = config.transition_length,
        Te_L = config.anode_Tev,
        Te_R = config.cathode_Tev,
        implicit_energy = config.implicit_energy,
        ion_temperature_K = config.ion_temperature_K,
        neutral_velocity = config.neutral_velocity,
        anode_mass_flow_rate = config.anode_mass_flow_rate,
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
