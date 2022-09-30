Base.@kwdef struct ElasticCollision{I} <: Reaction
    species::Species
    rate_coeff::I
end

abstract type ElectronNeutralModel <: ReactionModel end

struct NoElectronNeutral <: ElectronNeutralModel end

supported_species(::NoElectronNeutral) = Gas[]

load_reactions(::NoElectronNeutral, species) = ElasticCollision{Nothing}[]


struct LandmarkElectronNeutral <: ElectronNeutralModel end

supported_species(::LandmarkElectronNeutral) = [Xenon]

function load_reactions(::LandmarkElectronNeutral, species)
    landmark_electron_neutral = Returns(2.5e-13)
    return [ElasticCollision(Xenon(0), landmark_electron_neutral)]
end


struct GKElectronNeutral <: ElectronNeutralModel end

supported_species(::GKElectronNeutral) = [Xenon]

"""
    σ_en(Tev)

Electron neutral collision cross section in m² as a function of electron temperature in eV.
Eq. 3.6-13, from Fundamentals of Electric Propulsion, Goebel and Katz, 2008.

"""
@inline function σ_en(Tev)
    return max(0.0, 6.6e-19 * ((Tev / 4 - 0.1) / (1 + (Tev / 4)^1.6)))
end

function load_reactions(::GKElectronNeutral, species)
    rate_coeff(ϵ) = σ_en(2/3 * ϵ) * electron_sound_speed(2/3 * ϵ)
    return [ElasticCollision(Xenon(0), rate_coeff)]
end


Base.@kwdef struct ElectronNeutralLookup <: ElectronNeutralModel
    directories::Vector{String} = String[]
end

supported_species(::ElectronNeutralLookup) = Gas[]

function load_reactions(model::ElectronNeutralLookup, species)
    species_sorted = sort(species; by=x -> x.Z)
    reactions = ElasticCollision{LinearInterpolation{LinRangeWrapper{Float64, Int},Vector{Float64}}}[]
    folders = [model.directories; REACTION_FOLDER]
    product = nothing
    collision_type = "elastic"
    for i in 1:length(species)
        reactant = species_sorted[i]
        if reactant.Z > 0
            break
        else
            for folder in folders
                filename = rate_coeff_filename(reactant, product, collision_type, folder)
                if ispath(filename)
                    _, rate_coeff = load_rate_coeffs(reactant, product, collision_type,folder)
                    reaction = HallThruster.ElasticCollision(reactant, rate_coeff)
                    push!(reactions, reaction)
                    break
                end
            end
        end
    end

    return reactions
end
