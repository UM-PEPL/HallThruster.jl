struct ExcitationReaction <: Reaction
    reactant::Species
    rate_coeffs::Vector{Float64}
    energy::Float64
    function ExcitationReaction(energy::Float64, reactant::Species, rate_coeffs::Vector{Float64})
        return new(reactant, rate_coeffs, energy)
    end
end

LandmarkExcitationLookup() = :Landmark
ExcitationLookup() = :Lookup

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

function load_excitation_reactions(
        model::Symbol, species; directories = String[], kwargs...,
    )
    if model == :None
        return ExcitationReaction[]
    elseif model == :Landmark
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
    elseif model == :OVS
        # Rate coefficient for order verification studies
        Es = 0:1.0:255
        ks = ovs_rate_coeff_ex.(Es)
        return [ExcitationReaction(8.32, Xenon(0), ks)]
    elseif model == :Lookup
        species_sorted = sort(species; by = x -> x.Z)
        reactions = ExcitationReaction[]
        folders = [directories; REACTION_FOLDER]
        product = nothing
        collision_type = "excitation"

        for i in 1:length(species)
            reactant = species_sorted[i]

            filename = rate_coeff_filename(reactant, product, collision_type, nothing)
            filepath = find_file_in_dirs(filename, folders, cwd = false)

            if isnothing(filepath)
                continue
            end

            energy, rate_coeff = load_rate_coeff_file(filepath, collision_type)
            reaction = HallThruster.ExcitationReaction(energy, reactant, rate_coeff)
            push!(reactions, reaction)
        end
    else
        throw(ArgumentError("Invalid excitation model $(model). Choose :None, :Landmark, or :Lookup."))
    end
    return reactions
end
