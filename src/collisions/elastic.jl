Base.@kwdef struct ElasticCollision <: Reaction
    species::Species
    rate_coeffs::Vector{Float64}
end

NoElectronNeutral() = :None
ElectronNeutralLookup() = :Lookup
LandmarkElectronNeutral() = :Landmark

function load_elastic_collisions(model::Symbol, species; directories = String[], kwargs...)
    if model == :None
        return ElasticCollision[]
    elseif model == :Landmark
        return [ElasticCollision(Xenon(0), [2.5e-13, 2.5e-13, 2.5e-13])]
    elseif model == :Lookup
        species_sorted = sort(species; by = x -> x.Z)
        reactions = ElasticCollision[]
        folders = [directories; REACTION_FOLDER]
        product = nothing
        collision_type = "elastic"
        for i in 1:length(species)
            reactant = species_sorted[i]
            if reactant.Z > 0
                break
            end
            found = false
            for folder in folders
                filename = rate_coeff_filename(reactant, product, collision_type, folder)
                if ispath(filename)
                    _, rate_coeff = load_rate_coeffs(
                        reactant, product, collision_type, folder,
                    )
                    reaction = HallThruster.ElasticCollision(reactant, rate_coeff)
                    push!(reactions, reaction)
                    found = true
                    break
                end
            end
        end
        return reactions
    else
        throw(ArgumentError("Invalid elastic collision model $(model). Choose :None, :Landmark, or :Lookup."))
    end
end
