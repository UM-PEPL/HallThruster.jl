"""
    Gas
A chemical element in the gaseous state. Container for element properties used in fluid computations.

# Fields
`name::String`        Full name of gas (i.e. Xenon)

`short_name::String`  Short name/symbol (i.e. Xe for Xenon)

`γ::Float64`          Specific heat ratio / adiabatic index

`M::Float64`          Molar mass (grams/mol) or atomic mass units

`m::Float64`          Mass of atom in kg

`cp::Float64`         Specific heat at constant pressure in J / kg / K

`cv::Float64`         Specific heat at constant volume in J / kg / K

`R::Float64`          Gas constant in J / kg / K

"""
struct Gas
    name::String
    short_name::String
    γ::Float64      # Specific heat ratio / adiabatic index
    M::Float64      # Molar mass (grams/mol) or atomic mass units
    m::Float64      # Mass of atom in kg
    cp::Float64     # Specific heat at constant pressure
    cv::Float64     # Specific heat at constant volume
    R::Float64      # Gas constant
end

Base.show(io::IO, g::Gas) = print(io, g.name)
Base.show(io::IO, m::MIME"text/plain", g::Gas) = show(io, g)

"""
    Gas(name::String, short_name::String; γ::Float64, M::Float64)
Instantiate a new Gas, providing a name, short name, the adiabatic index, and the molar mass.
Other gas properties, including gas constant, specific heats at constant pressure/volume, and
mass of atom/molecule in kg will are then computed.

```jldoctest;setup = :(using HallThruster: Gas)
julia> Gas("Xenon", "Xe", γ = 5/3, M = 83.798)
Xenon
```
"""
function Gas(name, short_name; γ, M)
    R = R0 / M
    m = M / NA
    cp = γ / (γ - 1) * R
    cv = cp - R

    return Gas(name, short_name, γ, M, m, cp, cv, R)
end

# lets you do things like Xenon(1) == Species(Xenon, 1)
(g::Gas)(Z::Int) = Species(g, Z)

"""
    Species
Represents a gas with a specific charge state. In a plasma, different ionization states of the same
gas may coexist, so we need to be able to differentiate between these.

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
    element::Gas
    Z::Int
    symbol::Symbol
end

function Species(element::Gas, Z::Int)
    return Species(element, Z, Symbol(species_string(element, Z)))
end

Base.show(io::IO, s::Species) = print(io, string(s))
Base.show(io::IO, m::MIME"text/plain", s::Species) = show(io, s)

function species_string(element::Gas, Z::Int)
    sign_str = Z > 0 ? "+" : Z < 0 ? "-" : ""
    sign_str = abs(Z) > 1 ? "$(Z)" * sign_str : sign_str
    return element.short_name * sign_str
end

Base.string(s::Species) = string(s.symbol)

"""
    Air::Gas
Earth air at standard temperature and pressure
"""
const Air = Gas("Air", "Air"; γ=1.4, M=28.97)

"""
    Argon::Gas
Argon gas
"""
const Argon = Gas("Argon", "Ar"; γ=5 / 3, M=39.948)

"""
    Krypton::Gas
Krypton gas
"""
const Krypton = Gas("Krypton", "Kr"; γ=5 / 3, M=83.798)

"""
    Xenon::Gas
Xenon gas
"""
const Xenon = Gas("Xenon", "Xe"; γ=5 / 3, M=131.293)

"""
    Bismuth::Gas
Bismuth vapor
"""
const Bismuth = Gas("Bismuth", "Bi"; γ = 5/3, M = 208.9804)

"""
    Mercury::Gas
Mercury vapor
"""
const Mercury = Gas("Mercury", "Hg"; γ = 5/3, M = 200.59)

"""
    Electron::Species
Electron
"""
const Electron = Species(Gas("Electron", "e"; γ=5 / 3, M=1 / 1836.15267343), -1)
