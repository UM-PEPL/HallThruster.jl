
Base.@kwdef struct ConservationLawSystem
    type::Symbol
    nvars::Int
    u::Union{Float64,Nothing}
    T::Union{Float64,Nothing}
end

nvars(c::ConservationLawSystem) = c.nvars

function Base.show(io::IO, c::ConservationLawSystem)
    name = string(c.type)
    u = c.u
    T = c.T
    ustring = if isnothing(u)
        ""
    else
        "u = $u m/s"
    end
    Tstring = if isnothing(T)
        ""
    else
        "T = $T K"
    end
    argstring = join(filter(!isempty, [ustring, Tstring]), ", ")
    return print(io, name * "(" * argstring * ")")
end

"""
	ContinuityOnly
A `ConservationLawSystem` in which only continuity (mass conservation) is solved, while
velocity and temperature are held constant. Must specify a constant velocity (in m/s) and temperature (in K).

```jldoctest;setup = :(using HallThruster: ConservationLawSystem, ContinuityOnly)
julia> equation = ContinuityOnly(u = 300, T = 500)
ContinuityOnly(u = 300.0 m/s, T = 500.0 K)
```
"""
function ContinuityOnly(; u, T)
    return ConservationLawSystem(; type=:ContinuityOnly, nvars=1, u=Float64(u),
                                 T=Float64(T))
end
ContinuityOnly(u, T) = ContinuityOnly(; u, T)

"""
	IsothermalEuler
A `ConservationLawSystem` in which only continuity and inviscid momentum are solved, while
temperature is held constant. Must specify a constant temperature (in K).

```jldoctest;setup = :(using HallThruster: ConservationLawSystem, IsothermalEuler)
julia> equation = IsothermalEuler(T = 500)
IsothermalEuler(T = 500.0 K)
```
"""
function IsothermalEuler(; T)
    return ConservationLawSystem(; type=:IsothermalEuler, nvars=2, u=nothing, T=Float64(T))
end
IsothermalEuler(T) = IsothermalEuler(; T)

"""
	EulerEquations
A `ConservationLawSystem` for the inviscid Navier-Stokes equations, better known as the Euler equations.
Velocity and temperature are variable, so the values held in the ConservationLawSystem are set to zero
and subsequently unused.

```jldoctest;setup = :(using HallThruster: ConservationLawSystem, EulerEquations)
julia> equation = EulerEquations()
EulerEquations()
```
"""
function EulerEquations()
    return ConservationLawSystem(; type=:EulerEquations, nvars=3, u=nothing, T=nothing)
end