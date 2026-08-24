abstract type Reaction end

const EXCITATION_ENERGY_MERGE_TOLERANCE_EV = 1.0e-2

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
    fname = replace(
        fname, r"\((\d+)([+-]),(\d+)\*\)" => state -> begin
            parsed = match(r"\((\d+)([+-]),(\d+)\*\)", state)
            charge, sign, excited_level = parsed.captures
            return "$(charge)$(sign)_e$(excited_level)"
        end
    )

    # '*' is not legal in Windows filenames.
    fname = replace(
        fname, r"\((\d*)\*\)" => state -> begin
            parsed = match(r"\((\d*)\*\)", state)
            level = isempty(parsed.captures[1]) ? "1" : parsed.captures[1]
            return "_e$(level)"
        end
    )

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
            readline(io) # column header
        else
            energy = 0.0
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
    return rate_coeff(rxn, energy, ind)
end

function rate_coeff(rxn::Reaction, energy, ind::Int)
    isfinite(energy) || return first(rxn.rate_coeffs)
    N = length(rxn.rate_coeffs) - 2
    ind = ind > N ? N : ind < 0 ? 0 : ind
    r1 = rxn.rate_coeffs[ind + 1]
    r2 = rxn.rate_coeffs[ind + 2]
    return lerp(r1, r2, energy - ind)
end

"""Look up a rate using a previously clamped index and interpolation fraction."""
@inline function cached_rate_coeff(rxn::Reaction, ind::Int, fraction::Float64)
    r1 = rxn.rate_coeffs[ind + 1]
    r2 = rxn.rate_coeffs[ind + 2]
    return lerp(r1, r2, fraction)
end

"""Map each unique species symbol to its fluid-array index."""
function fluid_index_map(fluids)
    indices = Dict{Symbol, Int}()
    for (index, fluid) in enumerate(fluids)
        symbol = fluid.species.symbol
        if haskey(indices, symbol)
            throw(ArgumentError("Duplicate fluid species $symbol."))
        end
        indices[symbol] = index
    end
    return indices
end

function reaction_species_index(fluid_indices, species, role, reaction)
    index = get(fluid_indices, species.symbol, 0)
    index > 0 || throw(
        ArgumentError(
            "Missing $role species $(species.symbol) for reaction $reaction.",
        )
    )
    return index
end

reactant_indices(reactions, fluids::AbstractVector) =
    reactant_indices(reactions, fluid_index_map(fluids))

function reactant_indices(reactions, fluid_indices::AbstractDict)
    return [
        reaction_species_index(
                fluid_indices, reaction.reactant, "reactant", reaction,
            ) for reaction in reactions
    ]
end

product_indices(reactions, fluids::AbstractVector) =
    product_indices(reactions, fluid_index_map(fluids))

function product_indices(reactions, fluid_indices::AbstractDict)
    return [
        [
                reaction_species_index(fluid_indices, product, "product", reaction)
                for product in reaction.products
            ] for reaction in reactions
    ]
end

"""
    derive_species_energies(species, reactions)

Derive species energies in eV from one-to-one electron-impact reactions. The
neutral ground state defines zero energy for each gas. Reaction header energies
provide differences between levels, including ionization thresholds connecting
different charge-state manifolds.
"""
function derive_species_energies(species, reactions)
    adjacency = Dict{Symbol, Vector{Tuple{Symbol, Float64}}}(
        sp.symbol => Tuple{Symbol, Float64}[] for sp in species
    )

    for rxn in reactions
        length(rxn.products) == 1 || continue
        only(rxn.product_coeffs) == 1 || continue

        reactant = rxn.reactant
        product = only(rxn.products)
        reactant.element.formula == product.element.formula || continue

        push!(adjacency[reactant.symbol], (product.symbol, rxn.energy))
        push!(adjacency[product.symbol], (reactant.symbol, -rxn.energy))
    end

    energies = OrderedDict{Symbol, Float64}()
    pending = Symbol[]
    for sp in species
        if sp.Z == 0 && sp.excited_level == 0
            energies[sp.symbol] = 0.0
            push!(pending, sp.symbol)
        end
    end

    while !isempty(pending)
        source = popfirst!(pending)
        source_energy = energies[source]
        for (destination, delta_energy) in adjacency[source]
            candidate = source_energy + delta_energy
            if haskey(energies, destination)
                isapprox(
                    energies[destination], candidate;
                    rtol = 0.0, atol = EXCITATION_ENERGY_MERGE_TOLERANCE_EV,
                ) ||
                    error(
                    "Inconsistent excitation energy for $(destination): " *
                        "derived both $(energies[destination]) eV and $(candidate) eV " *
                        "from reaction headers (merge tolerance: " *
                        "$(EXCITATION_ENERGY_MERGE_TOLERANCE_EV) eV)."
                )
            else
                energies[destination] = candidate
                push!(pending, destination)
            end
        end
    end

    for sp in species
        if !haskey(energies, sp.symbol)
            error(
                "Could not derive the energy of species $(sp). Add a one-to-one " *
                    "excitation or ionization reaction connecting it to a species " *
                    "with a known energy."
            )
        end
    end

    return energies
end

function _configured_species(species_map, term, reaction)
    species_string = _species_string(term.species, term.charge, term.excited_level)
    species = get(species_map, Symbol(species_string), nothing)
    if isnothing(species) && term.excited_level == 0
        error("Species '$(species_string)' not found for reaction $(reaction).")
    end
    return species
end

function _configured_heavy_species(side, species_map, reaction)
    species = Species[]
    coefficients = UInt8[]
    for (term, coefficient) in side
        term.species == "e" && continue
        configured = _configured_species(species_map, term, reaction)
        # Reactions involving an excited level omitted from the propellant
        # configuration are omitted as a whole.
        isnothing(configured) && return nothing
        push!(species, configured)
        push!(coefficients, coefficient)
    end
    return species, coefficients
end

function _load_electron_impact_reaction(reaction, species_map, directories)
    lhs, rhs = _parse_reaction_equation(reaction["equation"])
    reactant_data = _configured_heavy_species(lhs, species_map, reaction)
    product_data = _configured_heavy_species(rhs, species_map, reaction)
    if isnothing(reactant_data) || isnothing(product_data)
        return nothing
    end
    reactants, reactant_coeffs = reactant_data
    products, product_coeffs = product_data

    length(reactants) == 1 || error(
        "Electron-impact reaction $(reaction) must have exactly one heavy-species reactant."
    )
    only(reactant_coeffs) == 1 || error(
        "Leading coefficient of species $(only(reactants)) must be one in reaction $(reaction)."
    )

    energy, rate_coeffs = _load_configured_rate_coefficients(
        reaction, "electron_impact", directories,
    )
    return ElectronImpactReaction(
        only(reactants), products, product_coeffs, rate_coeffs, energy,
    )
end

function _load_deexcitation_reaction(reaction, species_map)
    upper = _parse_species_term(reaction["target_species"])
    reactant = _configured_species(species_map, upper, reaction)
    isnothing(reactant) && return nothing

    levels = reaction["branches"]
    half_lives = reaction["half_lives"]
    length(levels) == length(half_lives) || error(
        "`branches` and `half_lives` must have equal length in de-excitation reaction $(reaction)."
    )

    products = Species[]
    rates = Float64[]
    for (level, half_life) in zip(levels, half_lives)
        lower = RxnTerm(upper.species, upper.charge, level)
        product = _configured_species(species_map, lower, reaction)
        isnothing(product) && continue
        push!(products, product)
        push!(rates, _deexcitation_rate(upper.excited_level, level, half_life, reaction))
    end
    isempty(products) && return nothing
    return DeExcitationReaction(reactant, products, rates)
end

function _load_target_collision(reaction, type, species_map, directories)
    target = _parse_species_term(reaction["target_species"])
    species = _configured_species(species_map, target, reaction)
    isnothing(species) && return nothing
    energy, rate_coeffs = _load_configured_rate_coefficients(reaction, type, directories)
    return type == "excitation" ?
        ExcitationReaction(energy, species, rate_coeffs) :
        ElasticCollision(species, rate_coeffs)
end

function load_reactions(
        propellant_config, species, iz_model, ex_model, en_model;
        directories = String[],
    )
    if !isempty(propellant_config) && isfile(propellant_config)
        contents = TOML.parsefile(propellant_config)
        if haskey(contents, "reactions")
            ei_reactions = ElectronImpactReaction[]
            ex_reactions = ExcitationReaction[]
            en_reactions = ElasticCollision[]
            de_reactions = DeExcitationReaction[]
            species_map = Dict{Symbol, Species}(s.symbol => s for s in species)

            for reaction in contents["reactions"]
                type = reaction["type"]
                if type == "electron_impact"
                    configured = _load_electron_impact_reaction(
                        reaction, species_map, directories,
                    )
                    isnothing(configured) || push!(ei_reactions, configured)
                elseif type == "de-excitation"
                    configured = _load_deexcitation_reaction(reaction, species_map)
                    isnothing(configured) || push!(de_reactions, configured)
                elseif type == "excitation" || type == "elastic"
                    configured = _load_target_collision(
                        reaction, type, species_map, directories,
                    )
                    if !isnothing(configured)
                        type == "excitation" ?
                            push!(ex_reactions, configured) :
                            push!(en_reactions, configured)
                    end
                else
                    error(
                        "Invalid reaction type $(type) in propellant config file " *
                            "$(propellant_config):\n$(reaction)"
                    )
                end
            end
            return ei_reactions, ex_reactions, en_reactions, de_reactions
        end
    end

    # Without configured reactions, generate the built-in collision set.
    return (
        load_electron_impact_reactions(iz_model, species; directories),
        load_excitation_reactions(ex_model, species; directories),
        load_elastic_collisions(en_model, species; directories),
        DeExcitationReaction[],
    )
end

function _load_configured_rate_coefficients(reaction, type, directories)
    rate_coeff_file = reaction["rate_coeff_file"]
    rate_coeff_path = find_file_in_dirs(rate_coeff_file, directories, cwd = true)
    isnothing(rate_coeff_path) && error(
        "Reaction rate coefficient file $(rate_coeff_file) not found in " *
            "provided directories $(directories)!"
    )
    return load_rate_coeff_file(rate_coeff_path, type)
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
