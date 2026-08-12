@public Gas, Species, Propellant

"""
	$(TYPEDEF)
A chemical element in the gaseous state. Container for element properties used in fluid computations.

# Fields
$(TYPEDFIELDS)
"""
struct Gas
    """Atomic symbol or molecular formula"""
    formula::Symbol
    """Specific heat ratio / adiabatic index"""
    γ::Float64
    """Molar mass (grams/mol) or atomic mass units"""
    M::Float64
    """Mass of atom in kg"""
    m::Float64
end

Base.show(io::IO, g::Gas) = print(io, string(g.formula))
Base.show(io::IO, ::MIME"text/plain", g::Gas) = show(io, g)

struct ElementTerm
    element::Symbol
    count::Int
end

struct MoleculeTerm
    components::Vector{Union{ElementTerm, MoleculeTerm}}
    count::Int
end

const ChemicalTerm = Union{ElementTerm, MoleculeTerm}

function parse_chemical_formula(formula::String)::Vector{ChemicalTerm}
    # a chemical formula has the following grammar (EBNF notation)
    # =============
    # nonzero_digit = "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9";
    # digit = "0" | nonzero_digit;
    # number = nonzero_digit, {digit};
    # element = [A-Z], [a-z];
    # element_with_number = element, [number]
    # molecule_with_parens = "(" + molecule + ")", [number]
    # molecule = element_with_number, {element_with_number | molecule_with_parens}

    parser = (;
        text = formula,
        pos = [1],
        token_start_pos = [1],
    )

    is_at_end(parser) = parser.pos[] > length(parser.text)

    peek(parser) = is_at_end(parser) ? '\0' : parser.text[parser.pos[]]

    function advance!(parser)
        parser.pos[] += 1
        return peek(parser)
    end

    function emit_token!(parser)
        tok = parser.text[parser.token_start_pos[]:(parser.pos[] - 1)]
        parser.token_start_pos[] = parser.pos[]
        return tok
    end

    function read_number!(parser)::Int
        while !is_at_end(parser) && isdigit(peek(parser))
            advance!(parser)
        end

        if parser.token_start_pos[] == parser.pos[]
            return 1
        else
            return parse(Int, emit_token!(parser))
        end
    end

    function expect!(parser, pred)
        @assert(pred(peek(parser)))
        return advance!(parser)
    end

    function consume_if!(parser, pred)
        if pred(peek(parser))
            return advance!(parser)
        else
            return peek(parser)
        end
    end

    function read_element!(parser)::ElementTerm
        # Consume first uppercase letter
        expect!(parser, isuppercase)

        # Read second lowercase letter if present
        consume_if!(parser, islowercase)
        elem_name = emit_token!(parser)

        # Next, attempt to read a digit. Return 1 if anything else found
        elem_count = read_number!(parser)
        return ElementTerm(Symbol(elem_name), elem_count)
    end

    function read_term!(parser)::ChemicalTerm
        if peek(parser) == '('
            advance!(parser)
            parser.token_start_pos[] = parser.pos[]
            molecule_components = read_molecule!(parser)
            if length(molecule_components) == 0
                throw(ErrorException("Parentheses in molecular formula '$(parser.text)' cannot be empty!"))
            end
            expect!(parser, ==(')'))
            parser.token_start_pos[] = parser.pos[]
            molecule_count = read_number!(parser)
            if length(molecule_components) == 1 && molecule_components[] isa ElementTerm
                return ElementTerm(
                    molecule_components[].element,
                    molecule_components[].count * molecule_count
                )
            else
                MoleculeTerm(molecule_components, molecule_count)
            end
        else
            return read_element!(parser)
        end
    end

    function read_molecule!(parser)::Vector{ChemicalTerm}
        terms = ChemicalTerm[]
        while !is_at_end(parser)
            c = peek(parser)
            if isuppercase(c) || c == '('
                term = read_term!(parser)
                push!(terms, term)
            else
                break
            end
        end
        return terms
    end

    return read_molecule!(parser)
end

function term_info(term)
    unit = if term isa ElementTerm
        (; mass = ELEMENTS[term.element].mass, num_atoms = 1)
    else
        molecule_info(term.components)
    end
    return (; mass = unit.mass * term.count, num_atoms = unit.num_atoms * term.count)
end

function molecule_info(components)
    num_atoms = 0
    mass = 0.0
    for term in components
        info = term_info(term)
        mass += info.mass
        num_atoms += info.num_atoms
    end
    return (; mass, num_atoms)
end

function term_formula(term)
    number_str(n) = n == 1 ? "" : "$(n)"

    if term isa ElementTerm
        return "$(string(term.element))" * number_str(term.count)
    else
        formula = molecule_formula(term.components)
        if term.count == 1
            return formula
        else
            return "(" * formula * ")" * number_str(term.count)
        end
    end
end

function molecule_formula(components)
    str = ""
    for term in components
        str *= term_formula(term)
    end
    return str
end

"""
    Gas(formula; γ, M)
Construct a gas from a chemical formula
For monatomic and diatomic species, the specific heat ratio is inferred
For triatomic and up, it must be provided
Unless provided, the full name will be determined from the short name provided.
"""
function Gas(formula::String; γ::Union{Float64, Nothing} = nothing, M::Union{Float64, Nothing} = nothing)
    components = parse_chemical_formula(formula)
    info = molecule_info(components)
    molecular_weight = M === nothing ? info.mass : M
    gamma = if γ === nothing
        if info.num_atoms == 1
            5 / 3
        elseif info.num_atoms == 2
            1.4
        else
            error("The ratio of specific heats, γ, must be provided for molecules with three or more atoms.")
        end
    else
        γ
    end

    formula = molecule_formula(components)
    return Gas(Symbol(formula), gamma, molecular_weight, molecular_weight / NA)
end

# lets you do things like Xenon(1) == Species(Xenon, 1)
(g::Gas)(Z::Int) = Species(g, Z)

"""
	$(TYPEDEF)
Represents a gas with a specific charge state. In a plasma, different ionization states of the same
gas may coexist, so we need to be able to differentiate between these.

# Fields
$(TYPEDFIELDS)

```jldoctest;setup = :(using HallThruster: Xenon, Species)
julia> Species(Xenon, 0)
Xe

julia> Species(Xenon, 1)
Xe(+)

julia> Species(Xenon, 3)
Xe(3+)
```
"""
struct Species
    """The gas that forms the base of the species"""
    element::Gas
    """The symbol of the species, i.e. `Symbol(Xe(+))` for `Species(Xenon, 1)`"""
    symbol::Symbol
    """The charge state of the species, i.e. Z = 1 for a singly-charged species"""
    Z::Int8
    @doc"""
    	Species(element::Gas, Z::Int) -> Species
    Construct a `Species` from a `Gas` and a charge state. You can also use the `(::Gas)(Z)` convenience constructor like so.

    ```julia
    julia> Xenon(0) == Species(Xenon, 0)
    true
    ```
    """ ->
    function Species(element::Gas, Z::Integer)::Species
        return new(element, Symbol(species_string(element, Z)), Int8(Z))
    end
end

Base.show(io::IO, s::Species) = print(io, string(s))
Base.show(io::IO, ::MIME"text/plain", s::Species) = show(io, s)

function _species_string(short_name::String, Z::Integer)
    sign_str = Z > 0 ? "+" : Z < 0 ? "-" : ""
    sign_str = abs(Z) > 1 ? "$(Z)" * sign_str : sign_str
    sign_str = Z > 0 ? "($(sign_str))" : ""
    return short_name * sign_str
end

species_string(element::Gas, Z::Integer) = _species_string(string(element.formula), Z)

Base.string(s::Species) = string(s.symbol)

# Built-in gases
# These are for convenience only and do not have any special status
const Hydrogen = Gas("H")
const MolecularHydrogen = Gas("H2")
const Helium = Gas("He")
const Oxygen = Gas("O")
const MolecularOxygen = Gas("O2")
const Nitrogen = Gas("N")
const MolecularNitrogen = Gas("N2")
const CarbonDioxide = Gas("CO2", γ = 1.28)
const Argon = Gas("Ar")
const Krypton = Gas("Kr")
const Xenon = Gas("Xe")
const Water = Gas("H2O", γ = 1.3)
const Bismuth = Gas("Bi")
const Mercury = Gas("Hg")
const Magnesium = Gas("Mg")


#=============================================================================
 Propellant
==============================================================================#
"""
    $(TYPEDEF)
Defines a propellant flowing through the thruster anode.
In addition to the neutral gas being used, the user specifies the anode mass flow rate and (optionally) the maximum charge state and temperature/velocity of the gas.

# Fields
$(TYPEDFIELDS)
"""
struct Propellant
    """
    A `Gas`. See [Propellants](propellants.md) for more. **Default:** `Xenon`.
    """
    gas::Gas
    """
    The mass flow rate of neutral atoms through the anode, in kg/s.
    """
    flow_rate_kg_s::Float64
    """
    Neutral velocity in m/s. **Default:** `$(DEFAULT_NEUTRAL_VELOCITY_M_S)`, or if `neutral_temperature` is set, that parameter is used to compute the velocity using a one-sided maxwellian flux approximation.
    """
    velocity_m_s::Float64
    """
    Neutral temperature in Kelvins for this propellant. **Default:** `$(DEFAULT_NEUTRAL_TEMPERATURE_K)`.
    """
    temperature_K::Float64
    """
    Ion temperature in Kelvins for this propellant. **Default:** `$(DEFAULT_ION_TEMPERATURE_K)`
    """
    ion_temperature_K::Float64
    """
    Allowed charges of ion states. **Default:** `[1, 2, ..., max_charge]`.
    """
    allowed_charges::Vector{Int}

    function Propellant(;
            gas::Union{String, Gas}, flow_rate_kg_s, max_charge = nothing,
            allowed_charges = nothing, velocity_m_s = nothing,
            temperature_K = nothing, ion_temperature_K = nothing,
        )

        if isnothing(velocity_m_s) && isnothing(temperature_K)
            # Use default values
            velocity_m_s = DEFAULT_NEUTRAL_VELOCITY_M_S
            temperature_K = DEFAULT_NEUTRAL_TEMPERATURE_K
        elseif isnothing(velocity_m_s)
            # Determine velocity from temperature
            temperature_K = convert_to_float64(temperature_K, units(:K))
            velocity_m_s = 0.25 * sqrt(8 * kB * temperature_K / π / gas.m)
        elseif isnothing(temperature_K)
            velocity_m_s = convert_to_float64(velocity_m_s, units(:m) / units(:s))
            temperature_K = DEFAULT_NEUTRAL_TEMPERATURE_K
        else
            velocity_m_s = convert_to_float64(velocity_m_s, units(:m) / units(:s))
            temperature_K = convert_to_float64(temperature_K, units(:K))
        end

        if isnothing(ion_temperature_K)
            ion_temperature_K = DEFAULT_ION_TEMPERATURE_K
        end

        final_allowed_charges =
        if !isnothing(max_charge) && !isnothing(allowed_charges)
            error("Provide either `max_charge` or `allowed_charges`, not both.")
        elseif !isnothing(allowed_charges)
            collect(allowed_charges)
        elseif !isnothing(max_charge)
            collect(1:max_charge)
        else
            [1]
        end

        length(unique(final_allowed_charges)) != length(final_allowed_charges) &&
            error("`allowed_charges` must not contain duplicates.")

        sort!(final_allowed_charges)

        ion_temperature_K = convert_to_float64(ion_temperature_K, units(:K))
        flow_rate_kg_s = convert_to_float64(flow_rate_kg_s, units(:kg) / units(:s))

        if gas isa String
            gas = Gas(gas)
        end

        return new(gas, flow_rate_kg_s, velocity_m_s, temperature_K, ion_temperature_K, final_allowed_charges)
    end
end

Propellant(gas; kwargs...) = Propellant(; gas, kwargs...)
Propellant(gas, flow_rate_kg_s; kwargs...) = Propellant(; gas, flow_rate_kg_s, kwargs...)

#=============================================================================
 Serialization
==============================================================================#
# Propellant list as of v0.21.7, the last version which saved propellants as strings alone
const propellants_v0_21_7 = Dict(
    "Xenon" => Xenon,
    "Krypton" => Krypton,
    "Argon" => Argon,
    "Bismuth" => Bismuth,
    "Mercury" => Mercury,
    "Nitrogen" => Nitrogen,
    "MolecularNitrogen" => MolecularNitrogen,
)

function Serialization.serialize(g::Gas)
    return OrderedDict("formula" => string(g.formula), "gamma" => g.γ)
end

function Serialization.deserialize(::Type{Gas}, d)
    if d isa String
        # backwards-compatible loading of gas from names alone for files written pre-v0.21.7
        return propellants_v0_21_7[d]
    else
        return Gas(d["formula"]; γ = d["gamma"])
    end
end

"""
    GASES
List of all built-in Gases. Users can of course provide their own if they want.
TODO: remove dependency on built-ins
"""
const GASES = [
    Argon, Krypton, Xenon, Nitrogen, MolecularNitrogen, Bismuth, Mercury, Water,
]
