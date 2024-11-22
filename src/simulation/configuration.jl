"""
$(TYPEDEF)
Hall thruster configuration struct. Only four mandatory fields: `discharge_voltage`, `thruster`, `anode_mass_flow_rate`, and `domain`.

# Fields
$(TYPEDFIELDS)
"""
struct Config{A <: AnomalousTransportModel, TC <: ThermalConductivityModel,
    W <: WallLossModel,
    S_N, S_IC, S_IM, S_ϕ, S_E,
    IC <: InitialCondition, HS <: HyperbolicScheme,}
    thruster::Thruster
    domain::Tuple{Float64, Float64}
    discharge_voltage::Float64
    anode_mass_flow_rate::Float64
    cathode_potential::Float64
    anode_Tev::Float64
    cathode_Tev::Float64
    wall_loss_model::W
    neutral_velocity::Float64
    neutral_temperature_K::Float64
    implicit_energy::Float64
    propellant::Gas
    ncharge::Int
    ion_temperature_K::Float64
    anom_model::A
    conductivity_model::TC
    ionization_model::Symbol
    excitation_model::Symbol
    electron_neutral_model::Symbol
    electron_ion_collisions::Bool
    min_number_density::Float64
    transition_length::Float64
    initial_condition::IC
    magnetic_field_scale::Float64
    source_neutrals::S_N
    source_ion_continuity::S_IC
    source_ion_momentum::S_IM
    source_potential::S_ϕ
    source_energy::S_E
    scheme::HS
    LANDMARK::Bool
    ion_wall_losses::Bool
    background_pressure_Torr::Float64
    background_temperature_K::Float64
    neutral_ingestion_multiplier::Float64
    anode_boundary_condition::Symbol
    anom_smoothing_iters::Int
    solve_plume::Bool
    apply_thrust_divergence_correction::Bool
    electron_plume_loss_scale::Float64
    reaction_rate_directories::Vector{String}
end

#=============================================================================
 Serialization of Config to JSON
==============================================================================#

# Don't write source terms to output or read them from input
function Serialization.exclude(::Type{C}) where {C <: Config}
    return (:source_neutrals, :source_ion_continuity,
        :source_ion_momentum, :source_potential, :source_energy,)
end

function Config(;
        thruster::Thruster,                 # MANDATORY ARGUMENT
        domain,                             # MANDATORY ARGUMENT
        discharge_voltage,                  # MANDATORY ARGUMENT
        anode_mass_flow_rate,               # MANDATORY ARGUMENT
        cathode_potential = 0.0,
        cathode_Tev = 3.0,
        anode_Tev = cathode_Tev,
        wall_loss_model::WallLossModel = WallSheath(BNSiO2, 1.0),
        neutral_velocity = nothing,
        neutral_temperature_K = nothing,
        implicit_energy::Number = 1.0,
        propellant::Gas = Xenon,
        ncharge::Int = 1,
        ion_temperature_K = 1000.0,
        anom_model::AnomalousTransportModel = TwoZoneBohm(1 / 160, 1 / 16),
        conductivity_model::ThermalConductivityModel = Mitchner(),
        ionization_model::Symbol = :Lookup,
        excitation_model::Symbol = :Lookup,
        electron_neutral_model::Symbol = :Lookup,
        electron_ion_collisions::Bool = true,
        min_number_density = 1e6,
        transition_length = 1e-3,
        initial_condition::IC = DefaultInitialization(),
        magnetic_field_scale::Float64 = 1.0,
        source_neutrals::S_N = nothing,
        source_ion_continuity::S_IC = nothing,
        source_ion_momentum::S_IM = nothing,
        source_potential::S_ϕ = Returns(0.0),
        source_energy::S_E = Returns(0.0),
        scheme::HyperbolicScheme = HyperbolicScheme(),
        LANDMARK = false,
        ion_wall_losses = false,
        background_pressure_Torr = 0.0,
        background_temperature_K = 100.0,
        neutral_ingestion_multiplier::Float64 = 1.0,
        anode_boundary_condition = :sheath,
        anom_smoothing_iters = 0,
        solve_plume = false,
        apply_thrust_divergence_correction = false,
        electron_plume_loss_scale = 1.0,
        reaction_rate_directories = String[],
) where {IC, S_N, S_IC, S_IM, S_ϕ, S_E}

    # check that number of ion source terms matches number of charges for both
    # continuity and momentum
    source_IC = ion_source_terms(ncharge, source_ion_continuity, "continuity")
    source_IM = ion_source_terms(ncharge, source_ion_momentum, "momentum")

    # Neutral source terms
    if isnothing(source_neutrals)
        source_neutrals = fill(Returns(0.0), 1)
    end

    # Convert to Float64 if using Unitful
    discharge_voltage = convert_to_float64(discharge_voltage, u"V")
    cathode_potential = convert_to_float64(cathode_potential, u"V")

    anode_Tev = convert_to_float64(anode_Tev, u"eV")
    cathode_Tev = convert_to_float64(cathode_Tev, u"eV")

    default_neutral_velocity = 150.0 # m/s
    default_neutral_temp = 500.0 # K
    if isnothing(neutral_velocity) && isnothing(neutral_temperature_K)
        neutral_velocity = default_neutral_velocity
        neutral_temperature_K = default_neutral_temp
    elseif isnothing(neutral_temperature_K)
        neutral_temperature_K = default_neutral_temp
        neutral_velocity = convert_to_float64(neutral_velocity, u"m/s")
    elseif isnothing(neutral_velocity)
        # compute neutral velocity from thermal speed
        neutral_temperature_K = convert_to_float64(neutral_temperature_K, u"K")
        neutral_velocity = 0.25 * sqrt(8 * kB * neutral_temperature_K / π / propellant.m)
    else
        neutral_velocity = convert_to_float64(neutral_velocity, u"m/s")
        neutral_temperature_K = convert_to_float64(neutral_temperature_K, u"K")
    end

    ion_temperature_K = convert_to_float64(ion_temperature_K, u"K")
    domain = (
        convert_to_float64(domain[1], u"m"),
        convert_to_float64(domain[2], u"m"),
    )

    anode_mass_flow_rate = convert_to_float64(anode_mass_flow_rate, u"kg/s")
    min_number_density = convert_to_float64(min_number_density, u"m^-3")

    background_temperature_K = convert_to_float64(
        background_temperature_K, u"K",)
    background_pressure_Torr = convert_to_float64(background_pressure_Torr, u"Pa")

    transition_length = convert_to_float64(transition_length, u"m")

    if anode_boundary_condition ∉ [:sheath, :dirichlet, :neumann]
        throw(ArgumentError("Anode boundary condition must be one of :sheath, :dirichlet, or :neumann. Got: $(anode_boundary_condition)"))
    end

    return Config(
        thruster, domain, discharge_voltage, anode_mass_flow_rate,
        cathode_potential, anode_Tev, cathode_Tev, wall_loss_model,
        neutral_velocity, neutral_temperature_K, implicit_energy, propellant, ncharge,
        ion_temperature_K, anom_model, conductivity_model,
        ionization_model, excitation_model, electron_neutral_model, electron_ion_collisions,
        min_number_density, transition_length,
        initial_condition, magnetic_field_scale, source_neutrals,
        source_IC,
        source_IM,
        source_potential,
        source_energy,
        scheme,
        LANDMARK,
        ion_wall_losses,
        background_pressure_Torr,
        background_temperature_K,
        neutral_ingestion_multiplier,
        anode_boundary_condition,
        anom_smoothing_iters,
        solve_plume,
        apply_thrust_divergence_correction,
        electron_plume_loss_scale,
        reaction_rate_directories,
    )
end

convert_to_float64(number::Number, unit) = Float64(number)
convert_to_float64(quantity::Quantity, unit) = uconvert(unit, quantity) |> ustrip |> Float64

function ion_source_terms(ncharge, source, type)
    if ncharge != length(source)
        throw(ArgumentError("Number of ion $type source terms must match number of charges"))
    end
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
        propellant(0); u = config.neutral_velocity, T = config.neutral_temperature_K,)
    ion_fluids = [IsothermalEuler(propellant(Z); T = config.ion_temperature_K)
                  for Z in 1:(config.ncharge)]

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
