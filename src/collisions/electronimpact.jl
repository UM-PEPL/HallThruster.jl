struct ElectronImpactReaction <: Reaction
    reactant::Species
    products::Vector{Species}
    product_coeffs::Vector{UInt8}
    rate_coeffs::Vector{Float64}
    energy::Float64
end

function ElectronImpactReaction(energy, reactant::Species, products::Vector{Species}, rate_coeffs)
    return ElectronImpactReaction(reactant, products, ones(UInt8, length(products)), rate_coeffs, energy)
end

LandmarkIonizationLookup() = :Landmark
IonizationLookup() = :Lookup

function Base.show(io::IO, rxn::ElectronImpactReaction)
    electron_input = "e(-)"
    net_charge = sum(prod.Z for prod in rxn.products) - rxn.reactant.Z
    electron_output = "$(net_charge == 0 ? "" : net_charge + 1)e(-)"
    reactant_str = string(rxn.reactant)
    rxn_str = electron_input * " + " * reactant_str * " -> "
    rxn_str *= electron_output
    for (prod, coeff) in zip(rxn.products, rxn.product_coeffs)
        product_str = string(prod)
        coeff_str = coeff == 1 ? "" : string(coeff)
        rxn_str *= " + $(coeff_str)$(product_str)"
    end
    return print(io, rxn_str)
end

function load_electron_impact_reactions(
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
        return [ElectronImpactReaction(12.12, Xenon(0), [Xenon(1)], rate_coeffs)]
    elseif model == :Lookup
        species_sorted = sort(species; by = x -> x.Z)
        reactions = ElectronImpactReaction[]
        folders = [directories; REACTION_FOLDER]

        # Check to make sure we find at least one reaction involving this species
        reactions_found = zeros(Bool, length(species))
        collision_type = "ionization"

        for i in 1:length(species)
            reactant = species_sorted[i]

            for j in (i + 1):length(species)
                product = species_sorted[j]

                filename = rate_coeff_filename(reactant, product, collision_type, nothing)
                filepath = find_file_in_dirs(filename, folders, cwd=false)

                if isnothing(filepath)
                    continue
                end

                energy, rate_coeff = load_rate_coeff_file(filepath, collision_type)
                reaction = ElectronImpactReaction(energy, reactant, [product], rate_coeff)

                push!(reactions, reaction)

                # Record that we found a reaction for both the reactant and product
                reactions_found[i] = true
                reactions_found[j] = true
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
        return [ElectronImpactReaction(12.12, Xenon(0), [Xenon(1)], [0.0, 0.0, 0.0])]
    else
        throw(ArgumentError("Invalid ionization model $(model). Select :Landmark or :Lookup"))
    end
end
