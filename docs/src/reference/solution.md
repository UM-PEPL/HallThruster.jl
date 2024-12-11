# Solution

```@meta
CurrentModule = HallThruster
```

## Types

```@docs
Solution
```

## Functions

```@docs
valid_fields()
alternate_field_names()
Base.getindex(sol::Solution, frame::Integer)
Base.getindex(sol::Solution, frames::AbstractVector)
Base.getindex(sol::Solution, field::Symbol)
Base.getindex(sol::Solution, field::Symbol, charge::Integer)
```
