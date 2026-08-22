# Outputs

```@meta
CurrentModule = HallThruster
```

Simulations in HallThruster.jl output `Solution` objects.
These contain the simulation's state at all requested timesteps in addition to the inputs the simulation was run with.
The outputs are organized into `Frame` objects, which contain arrays of data for all plasma properties of interest. Ground-state populations are grouped into neutral and ion `SpeciesState` objects, while explicitly tracked levels are available in `frame.excited_states` by full species symbol. Radiative transitions are stored in `frame.photon_emissions`, including their upper and lower states, spontaneous frequency, photon energy, and interval-averaged volumetric emission rate.

## Types

```@docs
Solution
Frame
SpeciesState
PhotonEmission
```

## Functions

```@docs
write_to_json
valid_fields
alternate_field_names
Base.getindex(sol::Solution, frame::Integer)
Base.getindex(sol::Solution, frames::AbstractVector)
```
