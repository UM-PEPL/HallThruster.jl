# Outputs

```@meta
CurrentModule = HallThruster
```

Simulations in HallThruster.jl output `Solution` objects.
These contain the simulation's state at all requested timesteps in addition to the inputs the simulation was run with.
The outputs are organized into `Frame` objects, which contain arrays of data for all plasma properties of interest, as well as `SpeciesState` objects, grouped by neutral species and ion species.

## Types

```@docs
Solution
Frame
SpeciesState
```

## Functions

```@docs
write_to_json
valid_fields
alternate_field_names
Base.getindex(sol::Solution, frame::Integer)
Base.getindex(sol::Solution, frames::AbstractVector)
```
