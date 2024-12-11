# Simulation parameters

```@meta
CurrentModule = HallThruster
```

In HallThruster, parameters used in the running of a specific simulation, such as the grid resolution and timestepping properties, are controlled by the `SimParams` struct.
Below are the fields of this struct.
As in all parts of the code, dimensional quantities are SI unless explicitly noted, but units may be provided using [`Unitful`](https://github.com/PainterQubits/Unitful.jl) or [`DynamicQuantities`](https://github.com/SymbolicML/DynamicQuantities.jl)

```@docs
SimParams{C <: CurrentController}
```

Using a `SimParams` in combination with a suitable [`Config`](../reference/config.md), we can use `run_simulation` to run a simulation.

```@docs
run_simulation(config::Config, sim::SimParams; kwargs...)
```
