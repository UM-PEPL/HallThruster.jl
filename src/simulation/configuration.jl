"""
$(TYPEDEF)
Hall thruster configuration struct. Only four mandatory fields: `discharge_voltage`, `thruster`, `anode_mass_flow_rate`, and `domain`.

# Fields
$(TYPEDFIELDS)
"""
@keywords mutable struct Config{
    A <: AnomalousTransportModel,
    TC <: ThermalConductivityModel,
    W <: WallLossModel,
    IC <: InitialCondition,
    HS <: HyperbolicScheme,
    # type parameters for source terms
    S_N,
    S_IC,
    S_IM,
    S_E,
}
    # Mandatory options
    thruster::Thruster
    discharge_voltage::Float64
    domain::Tuple{Float64, Float64}
    anode_mass_flow_rate::Float64

    # Optional options
    cathode_potential::Float64                = 0.0
    anode_Tev::Float64                        = 3.0
    cathode_Tev::Float64                      = 3.0
    anom_model::A                             = TwoZoneBohm(1 / 160, 1 / 16)
    conductivity_model::TC                    = Mitchner()
    wall_loss_model::W                        = WallSheath(BNSiO2, 1.0)
    neutral_velocity::Float64                 = 150.0
    neutral_temperature_K::Float64            = 500.0
    implicit_energy::Float64                  = 1.0
    propellant::Gas                           = Xenon
    ncharge::Int                              = 1
    ion_temperature_K::Float64                = 1000.0
    ionization_model::Symbol                  = :Lookup
    excitation_model::Symbol                  = :Lookup
    electron_neutral_model::Symbol            = :Lookup
    electron_ion_collisions::Bool             = true
    min_number_density::Float64               = 1e6
    transition_length::Float64                = 0.1
    initial_condition::IC                     = DefaultInitialization()
    scheme::HS                                = HyperbolicScheme()
    magnetic_field_scale::Float64             = 1.0
    source_neutrals::S_N                      = nothing
    source_ion_continuity::S_IC               = nothing
    source_ion_momentum::S_IM                 = nothing
    source_energy::S_E                        = nothing
    LANDMARK::Bool                            = false
    ion_wall_losses::Bool                     = false
    background_pressure_Torr::Float64         = 0.0
    background_neutral_temperature_K::Float64 = 150.0
    neutral_ingestion_multiplier::Float64     = 1.0
    anode_boundary_condition::Symbol          = :sheath
    anom_smoothing_iters::Int                 = 0
    solve_plume::Bool                         = false
    apply_thrust_divergence_correction::Bool  = false
    electron_plume_loss_scale::Float64        = 1.0
    reaction_rate_directories::Vector{String} = String[]
end

#=============================================================================
 Serialization of Config to JSON
==============================================================================#

StructTypes.StructType(::Config) = StructTypes.OrderedStruct()

# Don't write source terms to output
function StructTypes.excludes(::Type{C}) where {C <: Config}
    return (:source_neutrals, :source_ion_continuity,
        :source_ion_momentum, :source_potential, :source_energy,)
end

#=============================================================================#

function validate_config(config::Config)

    # check that anode BC selection is valid
    if config.anode_boundary_condition ∉ [:sheath, :dirichlet, :neumann]
        throw(ArgumentError(
            "Anode boundary condition must be one of :sheath, :dirichlet, or :neumann. \
            Got: $(config.anode_boundary_condition)",
        ))
    end

    # check that number of ion source terms matches number of charges for both
    # continuity and momentum
    @reset config.source_ion_continuity = check_ion_source_terms(
        config.ncharge, config.source_ion_continuity, "continuity",)
    @reset config.source_ion_momentum = check_ion_source_terms(
        config.ncharge, config.source_ion_momentum, "momentum",)

    # Set transition length default
    if config.transition_length === nothing
        @reset config.transition_length = 0.1 * config.thruster.geometry.channel_length
    end

    return config
end

function check_ion_source_terms(ncharge, source, type)
    if source === nothing
        return fill(nothing, ncharge)
    elseif ncharge != length(source)
        throw(
            ArgumentError("Number of ion $type source terms must match number of charges")
        )
    else
        return source
    end
end

function make_keys(fluid_range, subscript)
    len = length(fluid_range)
    if len == 1
        return (Symbol("ρ$(subscript)"))
    elseif len == 2
        return (Symbol("ρ$(subscript)"), Symbol("ρ$(subscript)u$(subscript)"))
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

function configure_fluids(config::Config)
    propellant = config.propellant

    neutral_fluid = ContinuityOnly(
        propellant(0); u = config.neutral_velocity, T = config.neutral_temperature_K,
    )
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
