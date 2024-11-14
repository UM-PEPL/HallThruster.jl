Base.@kwdef struct IonizationReaction <: Reaction
    energy::Float64
    reactant::Species
    product::Species
    rate_coeffs::Vector{Float64}
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

function load_ionization_reactions(model::Symbol, species; directories = String[], kwargs...)
    if model == :Lookup
        species_sorted = sort(species; by=x -> x.Z)
        reactions = IonizationReaction[]
        folders = [directories; REACTION_FOLDER]
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
                    throw(ArgumentError("No reactions including $(reactant) and $(product) in provided directories $(folders)."))
                end
            end
        end

        return reactions
    elseif model == :Landmark
        rates = readdlm(LANDMARK_RATES_FILE, ',', skipstart = 1)
        ϵ = rates[:, 1]
        k = rates[:, 2]
        itp = LinearInterpolation(ϵ, k)
        xs = 0:1.0:255
        rate_coeffs = itp.(xs)
        return [IonizationReaction(12.12, Xenon(0), Xenon(1), rate_coeffs)]
    else
        throw(ArgumentError("Invalid ionization moddl $(model)"))
    end
end
