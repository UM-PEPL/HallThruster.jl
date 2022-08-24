abstract type Reaction end

function rate_coeff_filename(reactant, product, reaction_type, folder = REACTION_FOLDER)
    fname = if product === nothing
        joinpath(folder, join([reaction_type, repr(reactant)], "_") * ".dat")
    else
        joinpath(folder, join([reaction_type, repr(reactant), repr(product)], "_") * ".dat")
    end
    return fname
end

@inline function maxwellian_vdf(Tev, v)
    sqrt(2/π) * v^2 * (me / e / Tev)^(3/2) * exp(-me * v^2 / 2 / e/ Tev)
end

function compute_rate_coeffs(energies, cross_section_func)
    integrand(Te, v) = let ϵ = 1/2 * me * v^2 / e
        cross_section_func(ϵ) * v * maxwellian_vdf(Te, v)
    end

    rate_coeffs = [
        quadgk(integrand $ (2/3 * ϵ), 0.0, 10 * sqrt(2 * e * ϵ / me))[1] for ϵ in energies
    ]
    return rate_coeffs
end

function compute_rate_coeffs(cross_section_filename)
    data = readdlm(cross_section_filename, ',')
    energies = 1.0:150.0

    ϵ, σ = data[:, 1], data[:, 2] * 1e-20
    cross_section_func = LinearInterpolation(ϵ, σ)
    rate_coeffs = compute_rate_coeffs(energies, cross_section_func)

    fname = splitpath(cross_section_filename)[end]
    open("reactions/$fname", "w") do io
        println(io, "Energy (eV)	Rate coefficient (m3/s)")
        writedlm(io, [energies rate_coeffs])
    end
end

function load_rate_coeffs(reactant, product, reaction_type, folder = REACTION_FOLDER)
    rates_file = rate_coeff_filename(reactant, product, reaction_type, folder)
    if ispath(rates_file)
        energy, rates = open(rates_file) do io
            if reaction_type == "ionization" || reaction_type == "excitation"
                energy = split(readline(io), ':')[2] |> strip |> parse $ (Float64)
            else
                energy = 0.0
            end
            rates = readdlm(io, skipstart=1)
            energy, rates
        end
    else
        throw(ArgumentError("$rates_file not found. This should normally be unreachable."))
    end
    ϵ = rates[:, 1]
    k = rates[:, 2]
    rate_coeff = LinearInterpolation(ϵ, k)
    return energy, rate_coeff
end

abstract type ReactionModel end

"""
    load_reactions(model::ReactionModel, species)::Vector{IonizationReaction}
Load ionization reactions for the provided species and ionization model
"""
@inline load_reactions(model::ReactionModel, species) = throw(ArgumentError("Function load_reactions($(typeof(model)), species) not implemented."))

function check_species(model::ReactionModel, species)
    supported = supported_gases(model)
    if length(supported) > 0
        for s in species
            if s.element ∉ supported
                throw(ArgumentError("$(s.element) is not supported by $(typeof(model)). The list of supported gases is $(supported)"))
            end
        end
    end
end

"""
    supported_gases(model::ReactionModel)::Vector{HallThruster.Gas}
Check which gases are supported by a given reaction model. If an empty vector is provided, then there are no restrictions on what gases can be used.
"""
@inline supported_gases(::ReactionModel) = Gas[]

"""
    maximum_charge_state(model::ReactionModel)::Int
Return the maximum supported charge state for a given reaction model. If 0 is returned, then no charge state restriction is applied.
"""
@inline maximum_charge_state(::ReactionModel) = 0

function check_charge_states(model::ReactionModel, species)
    max_charge = maximum_charge_state(model)
    if max_charge > 0
        for s in species
            if s.Z > max_charge
                throw(ArgumentError("The provided ionization model does not support ions with a charge state of $(s.Z). The maximum supported charge state is $(max_charge)."))
            end
        end
    end
end

function _load_reactions(model::ReactionModel, species)
    check_species(model, species)
    check_charge_states(model, species)
    load_reactions(model, species)
end

function _indices(symbol, reactions, species_range_dict)
    indices = [Int[] for _ in reactions]
    for (i, reaction) in enumerate(reactions)
        species = getfield(reaction, symbol).symbol
        ranges = species_range_dict[species]
        for r in ranges
            push!(indices[i], r[1])
        end
    end
    return indices
end

reactant_indices(reactions, species_range_dict) = _indices(:reactant, reactions, species_range_dict)
product_indices(reactions, species_range_dict) = _indices(:product, reactions, species_range_dict)
