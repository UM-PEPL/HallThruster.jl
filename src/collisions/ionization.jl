Base.@kwdef struct IonizationReaction <: Reaction
    energy::Float64
    reactant::Species
    product::Species
    rate_coeffs::Vector{Float64}
end

LandmarkIonizationLookup() = :Landmark
IonizationLookup() = :Lookup

function Base.show(io::IO, i::IonizationReaction)
    electron_input = "e-"
    electron_output = string(i.product.Z - i.reactant.Z + 1) * "e-"
    reactant_str = string(i.reactant)
    product_str = string(i.product)
    rxn_str = electron_input * " + " * reactant_str * " -> "
    rxn_str *= electron_output * " + " * product_str
    return print(io, rxn_str)
end

function load_ionization_reactions(
        model::Symbol, species; directories = String[], kwargs...,
    )
    if model == :Landmark
        # check species
        if length(species) > 2 || species[1] != Xenon(0) || species[2] != Xenon(1)
            throw(ArgumentError("Unsupported species $(species) for LANDMARK ionization lookup."))
        end

        rates = readdlm(LANDMARK_RATES_FILE, ',', skipstart = 1)
        ϵ = rates[:, 1]
        k = rates[:, 2]
        itp = LinearInterpolation(ϵ, k)
        xs = 0:1.0:255
        rate_coeffs = itp.(xs)
        return [IonizationReaction(12.12, Xenon(0), Xenon(1), rate_coeffs)]
    elseif model == :Lookup
        species_sorted = sort(species; by = x -> x.Z)
        reactions = IonizationReaction[]
        folders = [directories; REACTION_FOLDER]

        # Check to make sure we find at least one reaction involving this species
        reactions_found = zeros(Bool, length(species))

        for i in 1:length(species)
            reactant = species_sorted[i]

            for j in (i + 1):length(species)
                product = species_sorted[j]
                for folder in folders
                    filename = rate_coeff_filename(reactant, product, "ionization", folder)
                    if ispath(filename)
                        energy, rate_coeff = load_rate_coeffs(
                            reactant, product, "ionization", folder,
                        )
                        reaction = IonizationReaction(energy, reactant, product, rate_coeff)
                        push!(reactions, reaction)

                        # Record that we found a reaction for both the reactant and product
                        reactions_found[i] = true
                        reactions_found[j] = true
                        break
                    end
                end
            end
        end

        for i in eachindex(reactions_found)
            if !reactions_found[i]
                throw(ArgumentError("No reactions including $(species_sorted[i]) in provided directories $(folders)."))
            end
        end

        return reactions

    elseif model == :OVS
        # No ionization in OVS tests
        return [IonizationReaction(12.12, Xenon(0), Xenon(1), [0.0, 0.0, 0.0])]
    else
        throw(ArgumentError("Invalid ionization model $(model). Select :Landmark or :Lookup"))
    end
end
