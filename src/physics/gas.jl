@public Gas, Species

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
	short_name::String
	"""Specific heat ratio / adiabatic index"""
	γ::Float64
	"""Molar mass (grams/mol) or atomic mass units"""
	M::Float64
	"""Mass of atom in kg"""
	m::Float64
	"""Specific heat at constant pressure"""
	cp::Float64
	"""Specific heat at constant volume"""
	cv::Float64
	"""Gas constant"""
	R::Float64
	@doc"""
		Gas(name, short_name; γ, M) -> Gas
	Instantiate a new Gas, providing a name, short name, the adiabatic index, and the molar mass.
	Other gas properties, including gas constant, specific heats at constant pressure/volume, and
	mass of atom/molecule in kg will are then computed.

	```jldoctest;setup = :(using HallThruster: Gas)
	julia> Gas("Xenon", "Xe", γ = 5/3, M = 83.798)
	Xenon
	```
	"""->
	function Gas(name, short_name; γ, M)::Gas
		R = R0 / M
		m = M / NA
		cp = γ / (γ - 1) * R
		cv = cp - R

		return new(name, short_name, γ, M, m, cp, cv, R)
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
	"""The charge state of the species, i.e. Z = 1 for a singly-charged species"""
    Z::Int
	"""The symbol of the species, i.e. `Symbol(Xe+)` for `Species(Xenon, 1)`"""
    symbol::Symbol
	@doc"""
		Species(element::Gas, Z::Int) -> Species
	Construct a `Species` from a `Gas` and a charge state. You can also use the `(::Gas)(Z)` convenience constructor like so.

	```julia
	julia> Xenon(0) == Species(Xenon, 0)
	true
	```
	"""->
	function Species(element::Gas, Z::Int) :: Species
		return new(element, Z, Symbol(species_string(element, Z)))
	end
end

Base.show(io::IO, s::Species) = print(io, string(s))
Base.show(io::IO, ::MIME"text/plain", s::Species) = show(io, s)

function species_string(element::Gas, Z::Int)
    sign_str = Z > 0 ? "+" : Z < 0 ? "-" : ""
    sign_str = abs(Z) > 1 ? "$(Z)" * sign_str : sign_str
    return element.short_name * sign_str
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
	MolecularNitrogen::Gas
Molecular nitrogen gas
"""
const MolecularNitrogen = Gas("Molecular Nitrogen", "N2"; γ = 7/5, M = 28.0134)

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
 Serialization
==============================================================================#
const propellants = (; Xenon, Krypton, Argon, Bismuth, Mercury)
Serialization.SType(::Type{Gas}) = Serialization.Enum()
Serialization.options(::Type{Gas}) = propellants
