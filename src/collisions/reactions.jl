abstract type Reaction end


function rate_coeff_filename(reactant, product, reaction_type, folder = REACTION_FOLDER)
    fname = if product === nothing
        joinpath(folder, join([reaction_type, repr(reactant)], "_") * ".dat")
    else
        joinpath(folder, join([reaction_type, repr(reactant), repr(product)], "_") * ".dat")
    end
    return fname
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
    itp = LinearInterpolation(ϵ, k)
    xs = 0:1.0:255
    rate_coeffs = itp.(xs)
    return energy, rate_coeffs
end

abstract type ReactionModel end

"""
By default, rate_coeff looks for a lookup table stored in the reaction struct
"""
function rate_coeff(::ReactionModel, rxn::Reaction, energy)
    ind = Base.trunc(Int64, energy)
    N = length(rxn.rate_coeffs) - 2
    ind = ind > N ? N : ind
    r1 = rxn.rate_coeffs[ind+1]
    r2 = rxn.rate_coeffs[ind+2]
    t = energy - ind
    return (1 - t) * r1 + t * r2
end

"""
    load_reactions(model::ReactionModel, species)::Vector{IonizationReaction}
Load ionization reactions for the provided species and ionization model
"""
@inline load_reactions(model::ReactionModel, species; kwargs...) = throw(ArgumentError("Function load_reactions($(typeof(model)), species) not implemented."))

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

function _load_reactions(model::ReactionModel, species; kwargs...)
    check_species(model, species)
    check_charge_states(model, species)
    load_reactions(model, species; kwargs...)
end

function _indices(symbol, reactions, species_range_dict)
    indices = zeros(Int, length(reactions))
    for (i, reaction) in enumerate(reactions)
        species = getfield(reaction, symbol).symbol
        range = species_range_dict[species]
        indices[i] = range[1]
    end
    return indices
end

reactant_indices(reactions, species_range_dict) = _indices(:reactant, reactions, species_range_dict)
product_indices(reactions, species_range_dict) = _indices(:product, reactions, species_range_dict)
