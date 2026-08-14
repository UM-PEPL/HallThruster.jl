@public Gas, Species, Propellant

"""
	$(TYPEDEF)
A chemical element in the gaseous state. Container for element properties used in fluid computations.

# Fields
$(TYPEDFIELDS)
"""
struct Gas
    """Full name of gas (i.e. Xenon)"""
    name::String
    """Short name/symbol (i.e. Xe for Xenon)"""
    short_name::Symbol
    """Specific heat ratio / adiabatic index"""
    γ::Float64
    """Molar mass (grams/mol) or atomic mass units"""
    M::Float64
    """Mass of atom in kg"""
    m::Float64
    @doc"""
    	Gas(name, short_name; γ, M) -> Gas
    Instantiate a new Gas, providing a name, short name, the adiabatic index, and the molar mass.

    ```jldoctest;setup = :(using HallThruster: Gas)
    julia> Gas("Xenon", "Xe", γ = 5/3, M = 83.798)
    Xenon
    ```
    """ ->
    function Gas(name, short_name; γ, M)::Gas
        return new(name, Symbol(short_name), γ, M, M / NA)
    end
end

Base.show(io::IO, g::Gas) = print(io, g.name)
Base.show(io::IO, ::MIME"text/plain", g::Gas) = show(io, g)

# Convenience constructors such as Xenon(1) and Xenon(0, 1).
(g::Gas)(Z::Integer) = Species(g, Z)
(g::Gas)(Z::Integer, excited_level::Integer) = Species(g, Z, excited_level)

"""
	$(TYPEDEF)
Represents a gas with a specific charge state and electronic excitation level. In a plasma,
different ionization states of the same gas may coexist, so we need to be able to differentiate
between these. Excited states are likewise distinguished so they can be tracked as their own fluids.

# Fields
$(TYPEDFIELDS)

```jldoctest;setup = :(using HallThruster: Xenon, Species)
julia> Species(Xenon, 0)
Xe

julia> Species(Xenon, 1)
Xe(+)

julia> Species(Xenon, 3)
Xe(3+)

julia> Species(Xenon, 0, 1)
Xe(*)

julia> Species(Xenon, 0, 2)
Xe(2*)
```
"""
struct Species
    """The gas that forms the base of the species"""
    element::Gas
    """The symbol of the species, i.e. `Symbol(Xe(+))` for `Species(Xenon, 1)`"""
    symbol::Symbol
    """The charge state of the species, i.e. Z = 1 for a singly-charged species"""
    Z::Int8
    """The excitation level of the species; zero denotes the ground state"""
    excited_level::Int8
    @doc"""
        Species(element::Gas, Z::Integer, excited_level::Integer = 0) -> Species

    Construct a `Species` from a gas, charge state, and optional excitation level.
    An `excited_level` of zero represents the ground state.

    Call the gas directly with `gas(Z)` or `gas(Z, excited_level)` as a
    convenience:

    ```julia
    julia> Xenon(0) == Species(Xenon, 0)
    true

    julia> Xenon(0, 1) == Species(Xenon, 0, 1)
    true
    ```
    """ ->
    function Species(
            element::Gas,
            Z::Integer,
            excited_level::Integer = 0,
        )::Species
        excited_level < 0 && error("`excited_level` must be non-negative.")
        return new(
            element,
            Symbol(species_string(element, Z, excited_level)),
            Int8(Z),
            Int8(excited_level),
        )
    end
end

Base.show(io::IO, s::Species) = print(io, string(s))
Base.show(io::IO, ::MIME"text/plain", s::Species) = show(io, s)

"""
    is_excited(s::Species) -> Bool
Whether this species is electronically excited (`excited_level > 0`).
"""
is_excited(s::Species) = s.excited_level > 0

"""
    ground_state(s::Species) -> Species
The ground-state species with the same element and charge state as `s`.
"""
ground_state(s::Species) = Species(s.element, s.Z)

function _species_string(
        short_name::String, Z::Integer, excited_level::Integer = 0,
    )
    if Z == 0
        excited_level == 0 && return short_name
        level = excited_level == 1 ? "" : string(excited_level)
        return "$short_name($level*)"
    end

    sign = Z > 0 ? "+" : "-"
    if excited_level > 0
        return "$short_name($(abs(Z))$sign,$excited_level*)"
    end

    magnitude = abs(Z) == 1 ? "" : string(abs(Z))
    return "$short_name($magnitude$sign)"
end

species_string(element::Gas, Z::Integer, excited_level::Integer = 0) =
    _species_string(string(element.short_name), Z, excited_level)

Base.string(s::Species) = string(s.symbol)

"""
    Argon::Gas
Argon gas
"""
const Argon = Gas("Argon", "Ar"; γ = 5 / 3, M = 39.948)

"""
    Krypton::Gas
Krypton gas
"""
const Krypton = Gas("Krypton", "Kr"; γ = 5 / 3, M = 83.798)

"""
    Xenon::Gas
Xenon gas
"""
const Xenon = Gas("Xenon", "Xe"; γ = 5 / 3, M = 131.293)

"""
    Nitrogen::Gas
Atomic nitrogen gas
"""
const Nitrogen = Gas("Nitrogen", "N"; γ = 5 / 3, M = 14.0067)

"""
	MolecularNitrogen::Gas
Molecular nitrogen gas
"""
const MolecularNitrogen = Gas("Molecular Nitrogen", "N2"; γ = 7 / 5, M = 28.0134)

"""
    Bismuth::Gas
Bismuth vapor
"""
const Bismuth = Gas("Bismuth", "Bi"; γ = 5 / 3, M = 208.9804)

"""
    Mercury::Gas
Mercury vapor
"""
const Mercury = Gas("Mercury", "Hg"; γ = 5 / 3, M = 200.59)

"""
    GASES
List of all built-in Gases. Users can of course provide their own if they want.
"""
const GASES = [
    Argon, Krypton, Xenon, Nitrogen, MolecularNitrogen, Bismuth, Mercury,
]

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
    """
    Excitation levels tracked as their own neutral fluids (`Xe(*)`, `Xe(2*)`, ...), usable
    in reaction equations. **Default:** `[]` (excitation is lumped into an energy loss).
    """
    excited_levels::Vector{Int}
    """
    Excited levels tracked per ion charge state as their own fluids (`Xe(1+,1*)`, ...), keyed
    by charge. **Default:** empty (no excited ions).
    """
    excited_ion_levels::Dict{Int, Vector{Int}}

    function Propellant(;
            gas, flow_rate_kg_s, max_charge = nothing,
            allowed_charges = nothing, velocity_m_s = nothing,
            temperature_K = nothing, ion_temperature_K = nothing,
            excited_levels = nothing, excited_ion_levels = nothing,
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

        final_excited_levels = isnothing(excited_levels) ? Int[] : collect(excited_levels)

        any(<(1), final_excited_levels) &&
            error("`excited_levels` must be positive; the ground state (0) is always present.")

        length(unique(final_excited_levels)) != length(final_excited_levels) &&
            error("`excited_levels` must not contain duplicates.")

        sort!(final_excited_levels)

        final_excited_ion_levels = Dict{Int, Vector{Int}}()
        if !isnothing(excited_ion_levels)
            for (k, v) in excited_ion_levels
                Z = k isa AbstractString ? parse(Int, k) : Int(k)
                Z in final_allowed_charges ||
                    error("`excited_ion_levels` charge $(Z) is not in `allowed_charges`.")
                levels = collect(v)
                any(<(1), levels) &&
                    error("excited ion levels must be positive (ground state 0 is implicit).")
                length(unique(levels)) != length(levels) &&
                    error("excited ion levels must not contain duplicates.")
                final_excited_ion_levels[Z] = sort!(levels)
            end
        end

        ion_temperature_K = convert_to_float64(ion_temperature_K, units(:K))
        flow_rate_kg_s = convert_to_float64(flow_rate_kg_s, units(:kg) / units(:s))

        return new(
            gas, flow_rate_kg_s, velocity_m_s, temperature_K, ion_temperature_K,
            final_allowed_charges, final_excited_levels, final_excited_ion_levels,
        )
    end
end

Propellant(gas; kwargs...) = Propellant(; gas, kwargs...)
Propellant(gas, flow_rate_kg_s; kwargs...) = Propellant(; gas, flow_rate_kg_s, kwargs...)


#=============================================================================
 Serialization
==============================================================================#
const propellants = (; Xenon, Krypton, Argon, Bismuth, Mercury, MolecularNitrogen, Nitrogen)
Serialization.SType(::Type{Gas}) = Serialization.EnumType()
Serialization.options(::Type{Gas}) = propellants
