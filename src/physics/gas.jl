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
Xe+

julia> Species(Xenon, 3)
Xe3+
```
"""
struct Species
    """The gas that forms the base of the species"""
    element::Gas
    """The symbol of the species, i.e. `Symbol(Xe+)` for `Species(Xenon, 1)`"""
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

function species_string(element::Gas, Z::Integer)
    sign_str = Z > 0 ? "+" : Z < 0 ? "-" : ""
    sign_str = abs(Z) > 1 ? "$(Z)" * sign_str : sign_str
    return string(element.short_name) * sign_str
end

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
const Nitrogen = Gas("Nitrogen", "N"; γ = 5 / 3, M = 14.007)

"""
	MolecularNitrogen::Gas
Molecular nitrogen gas
"""
const MolecularNitrogen = Gas("Molecular Nitrogen", "N2"; γ = 7 / 5, M = 28.014)

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
    Maximum ion charge state. **Default:** 1.
    """
    max_charge::Int

    function Propellant(;
            gas, flow_rate_kg_s, max_charge = 1,
            velocity_m_s = nothing, temperature_K = nothing,
            ion_temperature_K = 1000.0,
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

        ion_temperature_K = convert_to_float64(ion_temperature_K, units(:K))
        flow_rate_kg_s = convert_to_float64(flow_rate_kg_s, units(:kg) / units(:s))

        return new(gas, flow_rate_kg_s, velocity_m_s, temperature_K, ion_temperature_K, max_charge)
    end
end

Propellant(gas; kwargs...) = Propellant(; gas, kwargs...)
Propellant(gas, flow_rate_kg_s; kwargs...) = Propellant(; gas, flow_rate_kg_s, kwargs...)

#=============================================================================
 Serialization
==============================================================================#
const propellants = (; Xenon, Krypton, Argon, Bismuth, Mercury)
Serialization.SType(::Type{Gas}) = Serialization.Enum()
Serialization.options(::Type{Gas}) = propellants
