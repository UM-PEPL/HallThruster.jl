Base.@kwdef struct ElasticCollision <: Reaction
    species::Species
    rate_coeffs::Vector{Float64}
end

abstract type ElectronNeutralModel <: ReactionModel end

# Definitions for NoElectronNeutral
struct NoElectronNeutral <: ElectronNeutralModel end
supported_species(::NoElectronNeutral) = Gas[]
load_reactions(::NoElectronNeutral, species) = ElasticCollision[]
rate_coeff(::NoElectronNeutral, ::Reaction, ::Float64) = 0.0

# Definitions for LandmarkElectronNeutral
struct LandmarkElectronNeutral <: ElectronNeutralModel end
supported_species(::LandmarkElectronNeutral) = [Xenon]
load_reactions(::LandmarkElectronNeutral, species) = [ElasticCollision(Xenon(0), Float64[])]
rate_coeff(::LandmarkElectronNeutral, ::Reaction, ::Float64) = 2.5e-13

# Definitions for GKElectronNeutral
"""
Electron neutral collision model as defined by
Eq. 3.6-13, from Fundamentals of Electric Propulsion, Goebel and Katz, 2008.
"""
struct GKElectronNeutral <: ElectronNeutralModel end
supported_species(::GKElectronNeutral) = [Xenon]
load_reactions(::GKElectronNeutral, species) = [ElasticCollision(Xenon(0), Float64[])]

@inline σ_en(Tev) = max(0.0, 6.6e-19 * ((Tev / 4 - 0.1) / (1 + (Tev / 4)^1.6)))

@inline function rate_coeff(::GKElectronNeutral, ::Reaction, energy::Float64)
    Tev = 2/3 * energy
    return σ_en(Tev) * electron_sound_speed(Tev)
end

# Definitions for ElectronNeutralLookup
Base.@kwdef struct ElectronNeutralLookup <: ElectronNeutralModel
    directories::Vector{String} = String[]
end

supported_species(::ElectronNeutralLookup) = Gas[]

function load_reactions(model::ElectronNeutralLookup, species)
    species_sorted = sort(species; by=x -> x.Z)
    reactions = []
    folders = [model.directories; REACTION_FOLDER]
    product = nothing
    collision_type = "elastic"
    for i in 1:length(species)
        reactant = species_sorted[i]
        if reactant.Z > 0
            break
        end
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

    return [r for r in reactions]
end
