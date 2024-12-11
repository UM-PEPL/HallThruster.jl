# Configuration

The `Config` struct contains all of the options you need specify the physics, geometry, and numerics of a simulation.
On this page, we will explain what options are available and what they do.
As in all parts of the code, dimensional quantities are SI unless explicitly noted, but units may be provided using [`Unitful`](https://github.com/PainterQubits/Unitful.jl) or [`DynamicQuantities`](https://github.com/SymbolicML/DynamicQuantities.jl)

```@meta
CurrentModule = HallThruster
```
```@docs
Config{A <: AnomalousTransportModel, TC <: ThermalConductivityModel, W <: WallLossModel, HS <: HyperbolicScheme, IC <: InitialCondition, S_N, S_IC, S_IM, S_E,}
```
