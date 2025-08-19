abstract type Reaction end

function rate_coeff_filename(reactant, product, reaction_type, folder = REACTION_FOLDER)
    fname = if product === nothing
        joinpath(folder, join([reaction_type, repr(reactant)], "_") * ".dat")
    else
        joinpath(folder, join([reaction_type, repr(reactant), repr(product)], "_") * ".dat")
    end
    return fname
end

function read_rate_coeff_file(path, reaction_type)
    if !ispath(path)
        throw(ArgumentError("Rate coefficient file $path not found."))
    end

    energy, rates = open(path) do io
        if reaction_type != "elastic"
            energy = parse(Float64, strip(split(readline(io), ':')[2]))
        else
            energy = 0.0
        end
        rates = readdlm(io, skipstart = 1)
        energy, rates
    end

    # Interpolate on grid from 0 to 255 eV of mean electron energy
    ϵ = rates[:, 1]
    k = rates[:, 2]
    itp = LinearInterpolation(ϵ, k)
    xs = 0:1.0:255
    rate_coeffs = itp.(xs)
    return energy, rate_coeffs
end

function load_rate_coeffs(reactant, product, reaction_type, folder = REACTION_FOLDER)
    rates_file = rate_coeff_filename(reactant, product, reaction_type, folder)
    return read_rate_coeff_file(rates_file, reaction_type)
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

function reactant_indices(reactions, fluids)
    indices = zeros(Int, length(reactions))
    for (i, reaction) in enumerate(reactions)
        species = reaction.reactant.symbol
        for (j, fluid) in enumerate(fluids)
            if fluid.species.symbol == species
                indices[i] = j
                break
            end
        end
    end
    return indices
end

function product_indices(reactions, fluids)
    indices = [Int[] for _ in eachindex(reactions)]
    for (i, reaction) in enumerate(reactions)
        for species in reaction.products
            for (j, fluid) in enumerate(fluids)
                if fluid.species.symbol == species.symbol
                    push!(indices[i], j)
                end
            end
        end
    end
    return indices
end
