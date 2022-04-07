Base.@kwdef struct IonizationReaction{I} <: Reaction
    energy::Float64
    reactant::Species
    product::Species
    rate_coeff::I
end

function Base.show(io::IO, i::IonizationReaction)
    electron_input = "e-"
    electron_output = string(i.product.Z - i.reactant.Z + 1) * "e-"
    reactant_str = string(i.reactant)
    product_str = string(i.product)
    rxn_str = electron_input * " + " * reactant_str * " -> "
    rxn_str *= electron_output * " + " * product_str
    return print(io, rxn_str)
end

abstract type IonizationModel <: ReactionModel end

"""
    IonizationLookup(;[directories::Vector{String} = String[]])
Default ionization model for HallThruster.jl.
Reads ionization rate coefficients from file. Looks (preferentially) in provided directories and in the reactions subfolder for rate coefficient files
"""
Base.@kwdef struct IonizationLookup <: IonizationModel
    directories::Vector{String} = String[]
end

@inline supported_gases(::IonizationLookup) = Gas[]
@inline maximum_charge_state(::IonizationLookup) = 0

function load_reactions(model::IonizationLookup, species)
    species_sorted = sort(species; by=x -> x.Z)
    reactions = IonizationReaction{LinearInterpolation{Float64,Float64}}[]
    folders = [model.directories; REACTION_FOLDER]
    for i in 1:length(species)
        for j in (i + 1):length(species)
            found = false
            reactant = species_sorted[i]
            product = species_sorted[j]
            for folder in folders
                filename = rate_coeff_filename(reactant, product, "ionization", folder)
                if ispath(filename)
                    energy, rate_coeff = load_rate_coeffs(reactant, product, "ionization",folder)
                    reaction = IonizationReaction(energy, reactant, product, rate_coeff)
                    push!(reactions, reaction)
                    found = true
                    break
                end
            end
            if !found
                throw(ArgumentError("No reactions including $(reactant) and $(product) in provided directories."))
            end
        end
    end

    return reactions
end


"""
    LandmarkIonizationLookup()
Ionization model for the LANDMARK benchmark.

Reads ionization rate coefficients from the landmark/landmark_rates.csv file in the HallThruster.jl main directory.
Supports only singly-charged Xenon.
"""
struct LandmarkIonizationLookup <: IonizationModel end

@inline supported_gases(::LandmarkIonizationLookup) = [Xenon]
@inline maximum_charge_state(::LandmarkIonizationLookup) = 1

function load_reactions(::LandmarkIonizationLookup, species)
    rates = readdlm(LANDMARK_RATES_FILE, ',', skipstart = 1)
    ϵ = rates[:, 1]
    k = rates[:, 2]
    rate_coeff = LinearInterpolation(ϵ, k)
    return [IonizationReaction(12.12, Xenon(0), Xenon(1), rate_coeff)]
end