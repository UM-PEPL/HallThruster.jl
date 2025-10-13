Base.@kwdef struct ElasticCollision <: Reaction
    reactant::Species
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

            filename = rate_coeff_filename(reactant, product, collision_type, nothing)
            filepath = find_file_in_dirs(filename, folders, cwd=false)

            if isnothing(filepath)
                continue
            end

            _, rate_coeff = load_rate_coeff_file(filepath, collision_type)
            reaction = HallThruster.ElasticCollision(reactant, rate_coeff)
            push!(reactions, reaction)
        end
        return reactions
    else
        throw(ArgumentError("Invalid elastic collision model $(model). Choose :None, :Landmark, or :Lookup."))
    end
end
