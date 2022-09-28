"""
$(TYPEDEF)
Hall thruster configuration struct. Only four mandatory fields: `discharge_voltage`, `thruster`, `anode_mass_flow_rate`, and `domain`.

# Fields
$(TYPEDFIELDS)
"""
struct Config{A<:AnomalousTransportModel, W<:WallLossModel, IZ<:IonizationModel, EX<:ExcitationModel, EN<:ElectronNeutralModel, HET<:Thruster, S_N, S_IC, S_IM, S_ϕ, S_E, T<:TransitionFunction, IC<:InitialCondition, CB, HS<:HyperbolicScheme}
    discharge_voltage::Float64
    cathode_potential::Float64
    anode_Te::Float64
    cathode_Te::Float64
    wall_loss_model::W
    neutral_velocity::Float64
    neutral_temperature::Float64
    implicit_energy::Float64
    propellant::Gas
    ncharge::Int
    ion_temperature::Float64
    anom_model::A
    ionization_model::IZ
    excitation_model::EX
    electron_neutral_model::EN
    electron_ion_collisions::Bool
    electron_pressure_coupled::Float64
    min_number_density::Float64
    min_electron_temperature::Float64
    transition_function::T
    initial_condition::IC
    callback::CB
    magnetic_field_scale::Float64
    source_neutrals::S_N
    source_ion_continuity::S_IC
    source_ion_momentum::S_IM
    source_potential::S_ϕ
    source_energy::S_E
    scheme::HS
    thruster::HET
    domain::Tuple{Float64, Float64}
    LANDMARK::Bool
    anode_mass_flow_rate::Float64
    ion_wall_losses::Bool
    solve_background_neutrals::Bool
    background_pressure::Float64
    background_neutral_temperature::Float64
    anode_boundary_condition::Symbol
end

function Config(;
        thruster::Thruster,                 # MANDATORY ARGUMENT
        domain,                             # MANDATORY ARGUMENT
        discharge_voltage,                  # MANDATORY ARGUMENT
        anode_mass_flow_rate,               # MANDATORY ARGUMENT
        cathode_potential                   = 0.0,
        cathode_Te                          = 3.0,
        anode_Te                            = cathode_Te,
        wall_loss_model::WallLossModel      = WallSheath(BNSiO2, 0.15),
        neutral_velocity                    = 300.0u"m/s",
        neutral_temperature                 = 300.0u"K",
        implicit_energy::Number             = 1.0,
        propellant::Gas                     = Xenon,
        ncharge::Int                        = 1,
        ion_temperature                     = 1000.0u"K",
        anom_model::AnomalousTransportModel = TwoZoneBohm(1/160, 1/16),
        ionization_model::IonizationModel   = IonizationLookup(),
        excitation_model::ExcitationModel   = ExcitationLookup(),
        electron_neutral_model::ElectronNeutralModel = ElectronNeutralLookup(),
        electron_ion_collisions::Bool       = true,
        electron_pressure_coupled::Number   = ncharge == 1 ? true : false,
        min_number_density                  = 1e6u"m^-3",
        min_electron_temperature            = min(anode_Te, cathode_Te),
        transition_function::TransitionFunction = LinearTransition(0.2 * thruster.geometry.channel_length, 0.0),
        initial_condition::IC               = DefaultInitialization(),
        callback                            = nothing,
        magnetic_field_scale::Float64       = 1.0,
        source_neutrals::S_N                = nothing,
        source_ion_continuity::S_IC         = nothing,
        source_ion_momentum::S_IM           = nothing,
        source_potential::S_ϕ               = Returns(0.0),
        source_energy::S_E                  = Returns(0.0),
        scheme::HyperbolicScheme            = HyperbolicScheme(),
        LANDMARK                            = false,
        ion_wall_losses                     = false,
        solve_background_neutrals           = false,
        background_pressure                 = 0.0u"Torr",
        background_neutral_temperature      = 100.0u"K",
        anode_boundary_condition            = :sheath
    ) where {IC, S_N, S_IC, S_IM, S_ϕ, S_E}

    # check that number of ion source terms matches number of charges for both
    # continuity and momentum
    source_IC = ion_source_terms(ncharge, source_ion_continuity, "continuity")
    source_IM = ion_source_terms(ncharge, source_ion_momentum,   "momentum")

    # Neutral source terms
    num_neutral_fluids = 1 + solve_background_neutrals
    if isnothing(source_neutrals)
        source_neutrals = fill(Returns(0.0), num_neutral_fluids)
    end

    # Convert to Float64 if using Unitful
    discharge_voltage = convert_to_float64(discharge_voltage, u"V")
    cathode_potential = convert_to_float64(cathode_potential, u"V")

    anode_Te = convert_to_float64(anode_Te, u"eV")
    cathode_Te = convert_to_float64(cathode_Te, u"eV")
    neutral_velocity = convert_to_float64(neutral_velocity, u"m/s")
    neutral_temperature = convert_to_float64(neutral_temperature, u"K")
    ion_temperature = convert_to_float64(ion_temperature, u"K")
    domain = (
        convert_to_float64(domain[1], u"m"),
        convert_to_float64(domain[2], u"m")
    )
    anode_mass_flow_rate = convert_to_float64(anode_mass_flow_rate, u"kg/s")
    min_electron_temperature = convert_to_float64(min_electron_temperature, u"eV")
    min_number_density = convert_to_float64(min_number_density, u"m^-3")

    background_neutral_temperature = convert_to_float64(background_neutral_temperature, u"K")
    background_pressure = convert_to_float64(background_pressure, u"Pa")

    if ncharge > 1 && electron_pressure_coupled > 0
        @warn("Electron pressure coupled method not compatible with multiply-charged ions. Switching to uncoupled scheme")
        electron_pressure_coupled = false
    end

    if anode_boundary_condition ∉ [:sheath, :dirichlet, :neumann]
        throw(ArgumentError("Anode boundary condition must be one of :sheath, :dirichlet, or :neumann. Got: $(anode_boundary_condition)"))
    end

    return Config(
        discharge_voltage, cathode_potential, anode_Te, cathode_Te, wall_loss_model,
        neutral_velocity, neutral_temperature, implicit_energy, propellant, ncharge, ion_temperature, anom_model,
        ionization_model, excitation_model, electron_neutral_model, electron_ion_collisions, Float64(electron_pressure_coupled), min_number_density, min_electron_temperature, transition_function,
        initial_condition, callback, magnetic_field_scale, source_neutrals,
        source_IC,
        source_IM,
        source_potential,
        source_energy,
        scheme,
        thruster,
        domain,
        LANDMARK,
        anode_mass_flow_rate,
        ion_wall_losses,
        solve_background_neutrals,
        background_pressure,
        background_neutral_temperature,
        anode_boundary_condition,
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
            Symbol("ρ$(subscript)u$(subscript)")
        )
    elseif len == 3
        return (
            Symbol("ρ$(subscript)"),
            Symbol("ρ$(subscript)u$(subscript)"),
            Symbol("ρ$(subscript)E$(subscript)")
        )
    else
        throw(ArgumentError("Too many equations on fluid (this should be unreachable)"))
    end
end

function configure_fluids(config)
    propellant = config.propellant

    anode_neutral_fluid = Fluid(propellant(0), ContinuityOnly(u = config.neutral_velocity, T = config.neutral_temperature))

    if config.solve_background_neutrals
        Tn_background = config.background_neutral_temperature
        background_neutral_fluid = Fluid(propellant(0), ContinuityOnly(u = background_neutral_velocity(config), T = Tn_background))
        neutral_fluids = [anode_neutral_fluid, background_neutral_fluid]
    else
        neutral_fluids = [anode_neutral_fluid]
    end

    ion_eqns = IsothermalEuler(T = config.ion_temperature)
    ion_fluids = [Fluid(propellant(Z), ion_eqns) for Z in 1:config.ncharge]

    fluids = [neutral_fluids; ion_fluids]

    species = [f.species for f in fluids]

    fluid_ranges = ranges(fluids)
    species_range_dict = Dict(Symbol(fluid.species) => UnitRange{Int64}[] for fluid in fluids)

    for (fluid, fluid_range) in zip(fluids, fluid_ranges)
        push!(species_range_dict[Symbol(fluid.species)], fluid_range)
    end

    return fluids, fluid_ranges, species, species_range_dict
end

function configure_index(fluids, fluid_ranges)
    lf = fluid_ranges[end][end]
    first_ion_fluid_index = findfirst(x -> x.species.Z > 0, fluids)

    keys_neutrals = (:ρn, )
    values_neutrals = (
        [f[1] for f in fluid_ranges[1:first_ion_fluid_index-1]],
    )

    keys_ions = (:ρi, :ρiui)
    values_ions = (
        [f[1] for f in fluid_ranges[first_ion_fluid_index:end]],
        [f[2] for f in fluid_ranges[first_ion_fluid_index:end]],
    )

    keys_fluids = (keys_neutrals..., keys_ions...)
    values_fluids = (values_neutrals..., values_ions...)
    keys_electrons = (:nϵ,)
    values_electrons = (lf + 1,)
    index_keys = (keys_fluids..., keys_electrons..., :lf)
    index_values = (values_fluids..., values_electrons..., lf)
    index = NamedTuple{index_keys}(index_values)
    return index
end
