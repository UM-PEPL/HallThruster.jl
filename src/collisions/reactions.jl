abstract type Reaction end

function rate_coeff_filename(reactant, product, reaction_type, folder = REACTION_FOLDER)
    fname = if product === nothing
        join([reaction_type, repr(reactant)], "_") * ".dat"
    else
        join([reaction_type, repr(reactant), repr(product)], "_") * ".dat"
    end

    # Remove '(' and ')' for backwards compatibility
    # TODO: switch all reaction file names to have the parentheses around charges for consistency
    fname = replace(fname, "(" => "", ")" => "")

    if !isnothing(folder)
        fname = joinpath(folder, fname)
    end

    return fname
end

function load_rate_coeff_file(path, reaction_type)
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

function load_reactions(propellant_config, species, iz_model, ex_model, en_model; directories = String[])
    if length(propellant_config) > 0 && isfile(propellant_config)
        contents = TOML.parsefile(propellant_config)

        if haskey(contents, "reactions")

            ei_reactions = ElectronImpactReaction[]
            ex_reactions = ExcitationReaction[]
            en_reactions = ElasticCollision[]

            species_map = Dict{Symbol, Species}(s.symbol => s for s in species)

            for reaction in contents["reactions"]

                type = reaction["type"]
                rate_coeff_file = reaction["rate_coeff_file"]
                rate_coeff_path = find_file_in_dirs(rate_coeff_file, directories, cwd=true)

                if isnothing(rate_coeff_path)
                    error("Reaction rate coefficient file $(rate_coeff_file) not found in provided directories $(directories)!")
                end

                energy, rate_coeffs = load_rate_coeff_file(rate_coeff_path, type)

                if type == "electron_impact"
                    lhs, rhs = _parse_reaction_equation(reaction["equation"])

                    reactants = Species[]
                    reactant_coeffs = UInt8[]
                    products = Species[]
                    product_coeffs = UInt8[]

                    for (side, species_arr, coeff_arr) in zip((lhs, rhs), (reactants, products), (reactant_coeffs, product_coeffs))
                        for (k, v) in side
                            if k.species == "e"
                                continue
                            end

                            target_species_str = _species_string(k.species, k.charge)
                            target_species_symbol = Symbol(target_species_str)
                            target_species = get(species_map, target_species_symbol, nothing)

                            if isnothing(target_species)
                                error("Species '$(target_species_str)' not found for reaction $(reaction).")
                            end

                            push!(species_arr, target_species)
                            push!(coeff_arr, v)
                        end 
                    end
                    
                    # Do some validation
                    if length(reactants) > 1
                        error("More than one reactant (excepting electrons) found for reaction $(reaction). Only single-reactant reactions are supported at present.")
                    end

                    if reactant_coeffs[1] != 1
                        error("Leading coefficient of species $(reactants[1]) must be one in reaction $(reaction).")
                    end

                    reaction = ElectronImpactReaction(reactants[1], products, product_coeffs, rate_coeffs, energy)
                    push!(ei_reactions, reaction)

                elseif type == "excitation" || type == "elastic"
                    target_species_str = reaction["target_species"]
                    target_species_symbol = Symbol(target_species_str)
                    target_species = get(species_map, target_species_symbol, nothing)
                    
                    if isnothing(target_species)
                        error("Species '$(target_species_str)' not found for reaction $(reaction).")
                    end

                    if type == "excitation"
                        push!(ex_reactions, ExcitationReaction(energy, target_species, rate_coeffs))
                    else
                        push!(en_reactions, ElasticCollision(target_species, rate_coeffs))
                    end
                else
                    error("Invalid reaction type $(type) in propellant config file $(propellant_config):\n$(reaction)")
                end
            end

            return ei_reactions, ex_reactions, en_reactions
        end
    end
    
    # If we're here, a file was not specified or there are no reactions in the file.
    ei_reactions = load_electron_impact_reactions(iz_model, species; directories)
    ex_reactions = load_excitation_reactions(ex_model, species; directories)
    en_reactions = load_elastic_collisions(en_model, species; directories)

    return ei_reactions, ex_reactions, en_reactions
end

#===========================================
 Reaction parsing utilities
============================================#

struct RxnTerm
    species::String
    charge::Int8
end

Base.show(io::IO, s::RxnTerm) = print(io, "$(_species_string(s.species, s.charge))")
Base.show(io::IO, ::MIME"text/plain", s::RxnTerm) = show(io, s)

mutable struct Lexer
    str::String
    index::Int
    function Lexer(s::String)
        return new(s, 1)
    end
end

function _peek(lex::Lexer)
    ci = nextind(lex.str, lex.index, 0)
    return lex.str[ci]
end

function _advance!(lex::Lexer)
    c = _peek(lex)
    lex.index = nextind(lex.str, lex.index, 1)
    return c
end

function _expect!(lex::Lexer, expected_char)
    got = _advance!(lex)
    if got != expected_char
        error("Expected $(expected_char) at position $(lex.index) in reaction equation $(lexer.str)")
    end
    return nothing
end

function _rest(lex::Lexer)
    return @view lex.str[lex.index:end]
end

function _takewhile!(pred, lex::Lexer)
    start_index = lex.index
    prev_index = 0
    for c in _rest(lex)
        if !pred(c)
            break
        end
        prev_index = lex.index
        lex.index += ncodeunits(c)
    end
    return @view lex.str[start_index:prev_index]
end

function _parse_number(lex)
    num_str = _takewhile!(isdigit, lex)
    if length(num_str) == 0
        return 1
    end
    return parse(Int, num_str)
end

function _is_term_char(char)
    return isdigit(char) || isletter(char)
end

function _parse_charge(lex)
    if lex.index > lastindex(lex.str) || _peek(lex) != '('
        return 0
    end

    _advance!(lex)
    charge = _parse_number(lex)

    sign = _advance!(lex)
    if sign == '-'
        charge = -charge
    elseif sign != '+'
        error("Expected + or - in sign term at position $(lex.index) in reaction equation $(lex.str)")
    end

    _expect!(lex, ')')

    return charge
end

function _parse_term!(lex)
    # Consume leading spaces
    _takewhile!(isspace, lex)

    # Get coefficient (set to one if not present)
    count = _parse_number(lex)

    # Get chemical symbol, erroring if one not found
    symbol = _takewhile!(_is_term_char, lex)

    if symbol == ""
        error("Expected chemical symbol in position $(lex.index) in reaction equation $(lex.str).")
    end

    # Get charge of species in parentheses. Set to zero if not present.
    charge = _parse_charge(lex)

    # Electron charge is -1 even if not specified
    if symbol == "e"
        charge = -1
    end

    # Consume trailing spaces
    _takewhile!(isspace, lex)

    return RxnTerm(symbol, charge), count
end

function _parse_side!(lex)
    terms = OrderedDict{RxnTerm, Int}()

    while lex.index < lastindex(lex.str) && _peek(lex) != '-'
        term, count = _parse_term!(lex)

        if haskey(terms, term)
            terms[term] += count
        else
            terms[term] = count
        end

        if lex.index >= lastindex(lex.str) || _peek(lex) == '-'
            break
        end

        _expect!(lex, '+')
    end

    return terms
end

function _parse_reaction_equation(eq::String)

    lex = Lexer(eq)
    lhs = _parse_side!(lex)
    _expect!(lex, '-')
    _expect!(lex, '>')
    rhs = _parse_side!(lex)

    # check total charge balance
    lhs_charge = sum(count * species.charge for (species, count) in pairs(lhs))
    rhs_charge = sum(count * species.charge for (species, count) in pairs(rhs))
    if lhs_charge != rhs_charge
        error("Charge does not balance in equation \"$(eq)\". Left charge: $(lhs_charge), right charge: $(rhs_charge).")
    end

    return lhs, rhs
end
