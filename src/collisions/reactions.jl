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
                energy = parse(Float64, strip(split(readline(io), ':')[2]))
            else
                energy = 0.0
            end
            rates = readdlm(io, skipstart = 1)
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

@inline lerp(a, b, t) = (1.0 - t) * a + t * b

"""
By default, rate_coeff looks for a lookup table stored in the reaction struct
"""
function rate_coeff(rxn::Reaction, energy)
    ind = Base.unsafe_trunc(Int, isfinite(energy) ? energy : 0)
    N = length(rxn.rate_coeffs) - 2
    ind = ind > N ? N : ind < 0 ? 0 : ind
    r1 = rxn.rate_coeffs[ind + 1]
    r2 = rxn.rate_coeffs[ind + 2]
    return lerp(r1, r2, energy - ind)
end

function _indices(symbol, reactions, fluids)
    indices = zeros(Int, length(reactions))

    for (i, reaction) in enumerate(reactions)
        species::Symbol = getfield(reaction, symbol).symbol
        for (j, fluid) in enumerate(fluids)
            if fluid.species.symbol == species
                indices[i] = j
                break
            end
        end
    end
    return indices
end

function reactant_indices(reactions, fluids)
    return _indices(:reactant, reactions, fluids)
end

function product_indices(reactions, fluids)
    return _indices(:product, reactions, fluids)
end
