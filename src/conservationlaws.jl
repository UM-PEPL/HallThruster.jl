abstract type ConservationLawSystem end

"""
	ContinuityOnly
A `ConservationLawSystem` in which only continuity (mass conservation) is solved, while
velocity and temperature are held constant. Must specify a constant velocity (in m/s) and temperature (in K).

```jldoctest;setup = :(using HallThruster: ContinuityOnly)
julia> equation = ContinuityOnly(u = 300, T = 500)
ContinuityOnly(300.0, 500.0)
```
"""
Base.@kwdef struct ContinuityOnly <: ConservationLawSystem
	u::Float64
	T::Float64
end

nvars(::Type{ContinuityOnly}) = 1

"""
	IsothermalEuler
A `ConservationLawSystem` in which only continuity and inviscid momentum are solved, while
temperature is held constant. Must specify a constant temperature (in K).

```jldoctest;setup = :(using HallThruster: IsothermalEuler)
julia> equation = IsothermalEuler(T = 500)
IsothermalEuler(500.0)
```
"""
Base.@kwdef struct IsothermalEuler <: ConservationLawSystem
	T::Float64
end

nvars(::Type{IsothermalEuler}) = 2

"""
	EulerEquations
A `ConservationLawSystem` for the inviscid Navier-Stokes equations, better known as the Euler equations

```jldoctest;setup = :(using HallThruster: EulerEquations)
julia> equation = EulerEquations()
EulerEquations()
```
"""
struct EulerEquations <: ConservationLawSystem end

nvars(::Type{EulerEquations}) = 3

