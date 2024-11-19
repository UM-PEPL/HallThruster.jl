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

# lets you do things like Xenon(1) == Species(Xenon, 1)
(g::Gas)(Z::Int) = Species(g, Z)

function Species(element::Gas, Z::Int)
    return Species(element, Z, Symbol(species_string(element, Z)))
end

Base.show(io::IO, s::Species) = print(io, string(s))
Base.show(io::IO, ::MIME"text/plain", s::Species) = show(io, s)

function species_string(element::Gas, Z::Int)
    sign_str = Z > 0 ? "+" : Z < 0 ? "-" : ""
    sign_str = abs(Z) > 1 ? "$(Z)" * sign_str : sign_str
    return element.short_name * sign_str
end

Base.string(s::Species) = string(s.symbol)
