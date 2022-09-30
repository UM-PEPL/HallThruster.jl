Base.@kwdef struct ExcitationReaction{I} <: Reaction
    energy::Float64
    reactant::Species
    rate_coeff::I
end

function Base.show(io::IO, i::ExcitationReaction)
    electron_input = "e-"
    electron_output = "e-"
    reactant_str = string(i.reactant)
    product_str = string(i.reactant) * "*"
    rxn_str = electron_input * " + " * reactant_str * " -> "
    rxn_str *= electron_output * " + " * product_str
    return print(io, rxn_str)
end

abstract type ExcitationModel <: ReactionModel end

"""
    NoExcitation <: ExcitationModel
Model for neglecting excitation energy losses
"""
struct NoExcitation <: ExcitationModel end

load_reactions(::NoExcitation, species) = ExcitationReaction{Nothing}[]

"""
    ExcitationLookup(;[directories::Vector{String} = String[]])
Default excitation model for HallThruster.jl.
Reads excitation rate coefficients from file. Looks (preferentially) in provided directories and in the reactions subfolder for rate coefficient files
"""
Base.@kwdef struct ExcitationLookup <: ExcitationModel
    directories::Vector{String} = String[]
end

function load_reactions(model::ExcitationLookup, species)
    species_sorted = sort(species; by=x -> x.Z)
    reactions = ExcitationReaction{LinearInterpolation{LinRangeWrapper{Float64, Int},Vector{Float64}}}[]
    folders = [model.directories; REACTION_FOLDER]
    product = nothing
    for i in 1:length(species)
        reactant = species_sorted[i]
        for folder in folders
            filename = rate_coeff_filename(reactant, product, "excitation", folder)
            if ispath(filename)
                energy, rate_coeff = load_rate_coeffs(reactant, product, "excitation",folder)
                reaction = HallThruster.ExcitationReaction(energy, reactant, rate_coeff)
                push!(reactions, reaction)
                break
            end
        end
    end

    return reactions
end


"""
    LandmarkExcitationLookup()
Excitation model for the LANDMARK benchmark.

Reads excitation rate coefficients from the landmark/landmark_rates.csv file in the HallThruster.jl main directory.
Supports only singly-charged Xenon.
"""
struct LandmarkExcitationLookup <: ExcitationModel end

function load_reactions(model::LandmarkExcitationLookup, species)
    rates = readdlm(LANDMARK_RATES_FILE, ',', skipstart = 1)
    ϵ = rates[:, 1]
    k_iz = rates[:, 2]
    k_loss = rates[:, 3]

    ionization_energy = 12.12
    excitation_energy = 8.32
    # We infer the excitation loss coefficient by subtracting the contribution of the ionization
    # loss coefficient from the inelastic loss term
    k_excitation = @. (k_loss - ionization_energy * k_iz) / excitation_energy
    rate_coeff = LinearInterpolation(ϵ, k_excitation)
    return [ExcitationReaction(excitation_energy, Xenon(0), rate_coeff)]
end
