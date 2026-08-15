abstract type Reaction end

struct DeExcitationReaction <: Reaction
    reactant::Species
    products::Vector{Species}
    rates::Vector{Float64}   # spontaneous emission rate per branch, 1/s
end

function rate_coeff_filename(reactant, product, reaction_type, folder = REACTION_FOLDER)
    fname = if product === nothing
        join([reaction_type, repr(reactant)], "_") * ".dat"
    else
        join([reaction_type, repr(reactant), repr(product)], "_") * ".dat"
    end

    # Charged excited states use Xe(2+,3*) -> Xe2+_e3.
    fname = replace(fname, r"\((\d+)([+-]),(\d+)\*\)" => state -> begin
        parsed = match(r"\((\d+)([+-]),(\d+)\*\)", state)
        charge, sign, excited_level = parsed.captures
        return "$(charge)$(sign)_e$(excited_level)"
    end)

    # '*' is not legal in Windows filenames.
    fname = replace(fname, r"\((\d*)\*\)" => state -> begin
        parsed = match(r"\((\d*)\*\)", state)
        level = isempty(parsed.captures[1]) ? "1" : parsed.captures[1]
        return "_e$(level)"
    end)

    if occursin('*', fname)
        error("Invalid excitation syntax in reaction filename: $(fname)")
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
        firstline = readline(io)
        if (reaction_type != "elastic") || (':' in firstline)
            energy = parse(Float64, strip(split(firstline, ':')[2]))
            header = readline(io)
        else
            energy = 0.0
            header = firstline
        end
        rates = readdlm(io)
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
    isfinite(energy) || return first(rxn.rate_coeffs)
    ind = Base.unsafe_trunc(Int, energy)
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
            de_reactions = DeExcitationReaction[]

            species_map = Dict{Symbol, Species}(s.symbol => s for s in species)

            for reaction in contents["reactions"]

                type = reaction["type"]

                # De-excitation is radiative: it has no rate coefficient file. `target_species`
                # is the upper state, `branches` the lower levels it decays to, and `half_lives`
                # the per-branch radiative half-lives.
                if type == "de-excitation"
                    upper_str = reaction["target_species"]
                    upper = _parse_species_term(upper_str)

                    reactant_str = _species_string(upper.species, upper.charge, upper.excited_level)
                    reactant = get(species_map, Symbol(reactant_str), nothing)
                    if isnothing(reactant)
                        error("Species '$(upper_str)' not found for de-excitation reaction $(reaction).")
                    end

                    levels = reaction["branches"]
                    half_lives = reaction["half_lives"]
                    if length(levels) != length(half_lives)
                        error("`branches` and `half_lives` must have equal length in de-excitation reaction $(reaction).")
                    end

                    products = Species[]
                    rates = Float64[]
                    for (level, half_life) in zip(levels, half_lives)
                        rate = _deexcitation_rate(
                            upper.excited_level, level, half_life, reaction,
                        )
                        product_str = _species_string(upper.species, upper.charge, level)
                        product = get(species_map, Symbol(product_str), nothing)
                        if isnothing(product)
                            setting = _excited_level_setting(upper.charge)
                            error(
                                "Product species '$(product_str)' not found for " *
                                    "de-excitation reaction $(reaction). Add level $(level) " *
                                    "to $(setting) on the corresponding propellant."
                            )
                        end
                        push!(products, product)
                        push!(rates, rate)
                    end

                    push!(de_reactions, DeExcitationReaction(reactant, products, rates))
                    continue
                end

                rate_coeff_file = reaction["rate_coeff_file"]
                rate_coeff_path = find_file_in_dirs(rate_coeff_file, directories, cwd = true)

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

                            target_species_str = _species_string(k.species, k.charge, k.excited_level)
                            target_species_symbol = Symbol(target_species_str)
                            target_species = get(species_map, target_species_symbol, nothing)

                            if isnothing(target_species)
                                if k.excited_level > 0
                                    setting = _excited_level_setting(k.charge)
                                    error(
                                        "Excited species '$(target_species_str)' not found for " *
                                            "reaction $(reaction). Add $(k.excited_level) to " *
                                            "$(setting) on the corresponding propellant."
                                    )
                                end
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
                    configured_target = reaction["target_species"]
                    target = _parse_species_term(configured_target)
                    target_species_str = _species_string(
                        target.species, target.charge, target.excited_level,
                    )
                    target_species = get(species_map, Symbol(target_species_str), nothing)

                    if isnothing(target_species)
                        error("Species '$(configured_target)' not found for reaction $(reaction).")
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

            return ei_reactions, ex_reactions, en_reactions, de_reactions
        end
    end

    # If we're here, a file was not specified or there are no reactions in the file.
    ei_reactions = load_electron_impact_reactions(iz_model, species; directories)
    ex_reactions = load_excitation_reactions(ex_model, species; directories)
    en_reactions = load_elastic_collisions(en_model, species; directories)
    de_reactions = DeExcitationReaction[]

    return ei_reactions, ex_reactions, en_reactions, de_reactions
end

#===========================================
 Reaction parsing utilities
============================================#

struct RxnTerm
    species::String
    charge::Int8
    excited_level::Int8
end

RxnTerm(species::String, charge) = RxnTerm(species, charge, 0)

Base.show(io::IO, s::RxnTerm) =
    print(io, _species_string(s.species, s.charge, s.excited_level))
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
        error("Expected $(expected_char) at position $(lex.index) in reaction equation $(lex.str)")
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

function _parse_excitation(lex)
    if lex.index > lastindex(lex.str) || _peek(lex) != '('
        return 0
    end

    suffix_start = lex.index
    _advance!(lex)
    level_str = _takewhile!(isdigit, lex)
    if lex.index > lastindex(lex.str) || _peek(lex) != '*'
        lex.index = suffix_start
        return 0
    end

    _advance!(lex)
    _expect!(lex, ')')
    excited_level = isempty(level_str) ? 1 : parse(Int, level_str)
    1 <= excited_level <= typemax(Int8) ||
        error("Invalid excitation level $(excited_level) in reaction equation $(lex.str)")
    return excited_level
end

function _parse_combined_state(lex)
    if lex.index > lastindex(lex.str) || _peek(lex) != '('
        return nothing
    end

    suffix_start = lex.index
    _advance!(lex)
    charge_str = _takewhile!(isdigit, lex)
    if isempty(charge_str)
        lex.index = suffix_start
        return nothing
    end
    charge_magnitude = parse(Int, charge_str)
    charge_magnitude > 0 ||
        error("Charge magnitude must be positive in reaction equation $(lex.str)")

    if lex.index > lastindex(lex.str) || !(_peek(lex) in ('+', '-'))
        lex.index = suffix_start
        return nothing
    end

    sign = _advance!(lex)
    charge = sign == '+' ? charge_magnitude : -charge_magnitude

    if lex.index > lastindex(lex.str) || _peek(lex) != ','
        lex.index = suffix_start
        return nothing
    end
    _advance!(lex)

    level_str = _takewhile!(isdigit, lex)
    isempty(level_str) &&
        error("Excitation level is required in reaction equation $(lex.str)")
    excited_level = parse(Int, level_str)
    if lex.index > lastindex(lex.str) || _peek(lex) != '*'
        lex.index = suffix_start
        return nothing
    end
    _advance!(lex)
    _expect!(lex, ')')

    if excited_level < 1 || excited_level > typemax(Int8)
        error(
            "Invalid excitation level $(excited_level) in reaction equation $(lex.str)"
        )
    end

    return charge, excited_level
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

    combined_state = _parse_combined_state(lex)
    if isnothing(combined_state)
        excited_level = _parse_excitation(lex)
        charge = _parse_charge(lex)
        excited_level > 0 && charge != 0 &&
            error("Use combined charged-excited notation in reaction equation $(lex.str)")
    else
        charge, excited_level = combined_state
    end

    # Electron charge is -1 even if not specified
    if symbol == "e"
        charge = -1
    end

    # Consume trailing spaces
    _takewhile!(isspace, lex)

    return RxnTerm(symbol, charge, excited_level), count
end

function _parse_species_term(text::AbstractString)
    lex = Lexer(String(text))
    term, count = _parse_term!(lex)
    count == 1 || error("Species term must not include a coefficient: $(text)")
    if lex.index <= lastindex(lex.str)
        error("Unexpected text at position $(lex.index) in species term $(text)")
    end
    return term
end

_excited_level_setting(charge) =
    charge == 0 ? "`excited_levels`" : "`excited_ion_levels[$charge]`"

function _deexcitation_rate(upper_level, lower_level, half_life, reaction)
    0 <= lower_level < upper_level ||
        error(
            "De-excitation branch $(upper_level) -> $(lower_level) must end at a " *
                "lower, non-negative excited level in reaction $(reaction)."
        )
    isfinite(half_life) && half_life > 0 ||
        error("De-excitation half-life must be positive and finite in reaction $(reaction).")
    return log(2.0) / half_life
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
