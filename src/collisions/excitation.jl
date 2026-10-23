@public ExcitationModel, ExcitationLookup, LandmarkExcitationLookup
@public NoExcitation, EffectiveProductionCost

"""
Abstract interface for models of electron energy losses due to excitation.
"""
abstract type ExcitationModel end

"""Load excitation reactions from rate-coefficient files."""
struct ExcitationLookup <: ExcitationModel end

"""Use excitation losses from the LANDMARK benchmark."""
struct LandmarkExcitationLookup <: ExcitationModel end

"""Disable excitation losses."""
struct NoExcitation <: ExcitationModel end

# Excitation model used by the order-verification tests.
struct OVSExcitation <: ExcitationModel end

"""
    EffectiveProductionCost([production_cost_multipliers]; A, B)

Use Fife's fit to Dugan's effective ion production cost model. Dictionary values
are the per-species ``C`` coefficients. `A` and `B` are also dictionaries keyed by
species. Omitted coefficients use Fife's Xenon values.
"""
struct EffectiveProductionCost <: ExcitationModel
    production_cost_multipliers::Dict{Symbol, Float64}
    A::Dict{Symbol, Float64}
    B::Dict{Symbol, Float64}
end

function production_cost_coefficient_dict(coefficients::AbstractDict)
    return Dict{Symbol, Float64}(
        (key isa Symbol ? key : Symbol(key)) => Float64(value)
            for (key, value) in coefficients
    )
end

function EffectiveProductionCost(
        multipliers::AbstractDict;
        A = Dict{Symbol, Float64}(), B = Dict{Symbol, Float64}(),
    )
    return EffectiveProductionCost(
        production_cost_coefficient_dict(multipliers),
        production_cost_coefficient_dict(A),
        production_cost_coefficient_dict(B),
    )
end

function EffectiveProductionCost(;
        production_cost_multipliers = Dict{Symbol, Float64}(),
        A = Dict{Symbol, Float64}(), B = Dict{Symbol, Float64}(),
    )
    return EffectiveProductionCost(production_cost_multipliers; A, B)
end

@inline function excitation_models()
    return (;
        ExcitationLookup,
        LandmarkExcitationLookup,
        NoExcitation,
        EffectiveProductionCost,
        OVSExcitation,
    )
end

function Serialization.SType(::Type{T}) where {T <: ExcitationModel}
    return Serialization.TaggedUnion()
end

Serialization.options(::Type{T}) where {T <: ExcitationModel} = excitation_models()

function excitation_model_from_symbol(model::Symbol)
    model == :Lookup && return ExcitationLookup()
    model == :Landmark && return LandmarkExcitationLookup()
    model == :None && return NoExcitation()
    model == :OVS && return OVSExcitation()
    model == :EffectiveProductionCost && return EffectiveProductionCost()
    throw(
        ArgumentError(
            "Invalid excitation model $(model). Choose :None, :Landmark, :Lookup, " *
                ":OVS, or :EffectiveProductionCost.",
        )
    )
end

normalize_excitation_model(model::ExcitationModel) = model
normalize_excitation_model(model::Symbol) = excitation_model_from_symbol(model)

is_effective_production_cost(::ExcitationModel) = false
is_effective_production_cost(::EffectiveProductionCost) = true

struct ExcitationRuntime
    use_effective_production_cost::Bool
    production_cost_multipliers::Dict{Symbol, Float64}
    A::Dict{Symbol, Float64}
    B::Dict{Symbol, Float64}
end

function excitation_runtime(::ExcitationModel)
    return ExcitationRuntime(
        false,
        Dict{Symbol, Float64}(),
        Dict{Symbol, Float64}(),
        Dict{Symbol, Float64}(),
    )
end

function excitation_runtime(model::EffectiveProductionCost)
    return ExcitationRuntime(
        true, model.production_cost_multipliers, model.A, model.B,
    )
end

struct ExcitationReaction <: Reaction
    reactant::Species
    rate_coeffs::Vector{Float64}
    energy::Float64
    function ExcitationReaction(energy::Float64, reactant::Species, rate_coeffs::Vector{Float64})
        return new(reactant, rate_coeffs, energy)
    end
end

const FIFE_XENON_PRODUCTION_COST_A = 0.254
const FIFE_XENON_PRODUCTION_COST_B = 0.677
const FIFE_XENON_PRODUCTION_COST_C = 2.0

@inline function effective_production_cost_multiplier(
        Tev, ionization_energy, A, B, C,
    )
    z = Tev / ionization_energy
    return A * exp(B / z) + C
end

function Base.show(io::IO, i::ExcitationReaction)
    electron_input = "e(-)"
    electron_output = "e(-)"
    reactant_str = string(i.reactant)
    product_str = string(i.reactant) * "(*)"
    rxn_str = electron_input * " + " * reactant_str * " -> "
    rxn_str *= electron_output * " + " * product_str
    return print(io, rxn_str)
end

ovs_rate_coeff_ex(ϵ) = 1.0e-12 * exp(-8.32 / ϵ)

load_excitation_reactions(::NoExcitation, _species; kwargs...) = ExcitationReaction[]
load_excitation_reactions(::EffectiveProductionCost, _species; kwargs...) =
    ExcitationReaction[]

function load_excitation_reactions(
        ::LandmarkExcitationLookup, _species; kwargs...,
    )
    rates = readdlm(LANDMARK_RATES_FILE, ',', skipstart = 1)
    ϵ = rates[:, 1]
    k_iz = rates[:, 2]
    k_loss = rates[:, 3]

    ionization_energy = 12.12
    excitation_energy = 8.32
    # We infer the excitation loss coefficient by subtracting the contribution of the ionization
    # loss coefficient from the inelastic loss term
    k_excitation = @. (k_loss - ionization_energy * k_iz) / excitation_energy
    itp = LinearInterpolation(ϵ, k_excitation)
    xs = 0:1.0:255
    rate_coeffs = itp.(xs)
    return [ExcitationReaction(excitation_energy, Xenon(0), rate_coeffs)]
end

function load_excitation_reactions(::OVSExcitation, _species; kwargs...)
    Es = 0:1.0:255
    ks = ovs_rate_coeff_ex.(Es)
    return [ExcitationReaction(8.32, Xenon(0), ks)]
end

function load_excitation_reactions(
        ::ExcitationLookup, species; directories = String[], kwargs...,
    )
    species_sorted = sort(species; by = x -> x.Z)
    reactions = ExcitationReaction[]
    folders = [directories; REACTION_FOLDER]
    product = nothing
    collision_type = "excitation"

    for reactant in species_sorted
        filename = rate_coeff_filename(reactant, product, collision_type, nothing)
        filepath = find_file_in_dirs(filename, folders, cwd = false)

        if isnothing(filepath)
            continue
        end

        energy, rate_coeff = load_rate_coeff_file(filepath, collision_type)
        reaction = HallThruster.ExcitationReaction(energy, reactant, rate_coeff)
        push!(reactions, reaction)
    end
    return reactions
end

function load_excitation_reactions(model::ExcitationModel, _species; kwargs...)
    throw(
        ArgumentError(
            "load_excitation_reactions is not implemented for $(typeof(model)).",
        )
    )
end

function load_excitation_reactions(model::Symbol, species; kwargs...)
    return load_excitation_reactions(excitation_model_from_symbol(model), species; kwargs...)
end
