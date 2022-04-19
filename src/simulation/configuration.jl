"""
$(TYPEDEF)
Hall thruster configuration struct. Only four mandatory fields: `discharge_voltage`, `thruster`, `anode_mass_flow_rate`, and `domain`.

# Fields
$(TYPEDFIELDS)
"""
struct Config{A<:AnomalousTransportModel, W<:WallLossModel, IZ<:IonizationModel, EX<:ExcitationModel, EN<:ElectronNeutralModel, HET<:Thruster, S_N, S_IC, S_IM, S_ϕ, S_E, T<:TransitionFunction, IC<:InitialCondition, CB, HS<:HyperbolicScheme, ECL}
    discharge_voltage::Float64
    cathode_potential::Float64
    anode_Te::Float64
    cathode_Te::Float64
    wall_loss_model::W
    wall_collision_freq::Float64
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
    progress_interval::Int
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
    electron_cond_lookup::ECL
end

function Config(;
        discharge_voltage,                  # MANDATORY ARGUMENT
        cathode_potential                   = 0.0,
        anode_Te                            = 3.0,
        cathode_Te                          = 3.0,
        wall_loss_model::WallLossModel      = ConstantSheathPotential(sheath_potential = -20.0, inner_loss_coeff = 1.0, outer_loss_coeff = 1.0),
        wall_collision_freq                 = 0.0,
        neutral_velocity                    = 300.0,
        neutral_temperature                 = 300.0,
        implicit_energy::Number             = 1.0,
        propellant::Gas                     = Xenon,
        ncharge::Int                        = 1,
        ion_temperature                     = 1000.,
        anom_model::AnomalousTransportModel = TwoZoneBohm(1/160, 1/16),
        ionization_model::IonizationModel   = IonizationLookup(),
        excitation_model::ExcitationModel   = ExcitationLookup(),
        electron_neutral_model::ElectronNeutralModel = ElectronNeutralLookup(),
        electron_ion_collisions::Bool       = true,
        electron_pressure_coupled::Number   = 1.0,
        min_number_density                  = 1e6,
        min_electron_temperature            = 1.0,
        transition_function::TransitionFunction = LinearTransition(0.001, 0.0),
        progress_interval::Int              = 0,
        initial_condition::IC               = DefaultInitialization(),
        callback                            = nothing,
        magnetic_field_scale::Float64       = 1.0,
        source_neutrals::S_N                = Returns(0.0),
        source_ion_continuity::S_IC         = nothing,
        source_ion_momentum::S_IM           = nothing,
        source_potential::S_ϕ               = Returns(0.0),
        source_energy::S_E                  = Returns(0.0),
        scheme::HyperbolicScheme            = HyperbolicScheme(flux_function = rusanov, limiter = minmod, reconstruct = false, WENO = false),
        thruster::Thruster,                 # MANDATORY ARGUMENT
        domain,                             # MANDATORY ARGUMENT
        LANDMARK                            = false,
        anode_mass_flow_rate,               # MANDATORY ARGUMENT,
        electron_cond_lookup                = ElectronCondLookup(),
    ) where {IC, S_N, S_IC, S_IM, S_ϕ, S_E}

    # check that number of ion source terms matches number of charges for both
    # continuity and momentum
    source_IC = ion_source_terms(ncharge, source_ion_continuity, "continuity")
    source_IM = ion_source_terms(ncharge, source_ion_momentum,   "momentum")

    # Convert to Float64 if using Unitful
    discharge_voltage = convert_to_float64(discharge_voltage, u"V")
    cathode_potential = convert_to_float64(cathode_potential, u"V")
    anode_Te = convert_to_float64(anode_Te, u"eV")
    cathode_Te = convert_to_float64(cathode_Te, u"eV")
    wall_collision_freq = convert_to_float64(wall_collision_freq, u"Hz")
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

    return Config(
        discharge_voltage, cathode_potential, anode_Te, cathode_Te, wall_loss_model, wall_collision_freq,
        neutral_velocity, neutral_temperature, implicit_energy, propellant, ncharge, ion_temperature, anom_model,
        ionization_model, excitation_model, electron_neutral_model, electron_ion_collisions, Float64(electron_pressure_coupled), min_number_density, min_electron_temperature, transition_function,
        progress_interval, initial_condition, callback, magnetic_field_scale, source_neutrals,
        source_IC, source_IM, source_potential, source_energy, scheme, thruster, domain, LANDMARK, anode_mass_flow_rate, electron_cond_lookup
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

function configure_index(fluid_ranges)
    lf = fluid_ranges[end][end]

    ncharge = length(fluid_ranges)-1
    solve_ion_temp = length(fluid_ranges[2]) == 3

    keys_neutrals = (:ρn, )
    values_neutrals = (1, )

    if solve_ion_temp
        keys_ions = (:ρi, :ρiui, :ρiuiEi)
        values_ions = (
            [f[1] for f in fluid_ranges[2:end]]...,
            [f[2] for f in fluid_ranges[2:end]]...,
            [f[3] for f in fluid_ranges[2:end]]...,
        )
    else
        keys_ions = (:ρi, :ρiui)
        values_ions = (
            [f[1] for f in fluid_ranges[2:end]],
            [f[2] for f in fluid_ranges[2:end]],
        )
    end

    keys_fluids = (keys_neutrals..., keys_ions...)
    values_fluids = (values_neutrals..., values_ions...)
    keys_electrons = (:nϵ, :Tev, :ne, :pe, :ϕ, :grad_ϕ, :ue)
    values_electrons = lf .+ collect(1:7)
    index_keys = (keys_fluids..., keys_electrons..., :lf)
    index_values = (values_fluids..., values_electrons..., lf)
    index = NamedTuple{index_keys}(index_values)
    return index
end

function configure_fluids(config)
    propellant = config.propellant
    species = [propellant(i) for i in 0:config.ncharge]
    neutral_fluid = Fluid(species[1], ContinuityOnly(u = config.neutral_velocity, T = config.neutral_temperature))
    ion_eqns = IsothermalEuler(T = config.ion_temperature)
    ion_fluids = [Fluid(species[i+1], ion_eqns) for i in 1:config.ncharge]
    fluids = [neutral_fluid; ion_fluids]
    fluid_ranges = ranges(fluids)
    species_range_dict = Dict(Symbol(fluid.species) => fluid_range
                              for (fluid, fluid_range) in zip(fluids, fluid_ranges))
    return fluids, fluid_ranges, species, species_range_dict
end