# Postprocessing

```@meta
CurrentModule = HallThruster
```
## Types

```@docs
Postprocess
```

## Functions

### Time-averaging

Time-averaging can be accomplished by providing a start time, or a starting frame.

```@docs
time_average(sol::Solution, start_time)
time_average(sol::Solution, start_frame::Integer = 1)
```

### Global metrics

These functions compute global metrics, i.e. thrust, currents, and efficiencies.
Each has a version that computes the metric for the entire solution, and a version that acts on a specific frame.

```@docs
thrust(sol::Solution, frame::Integer)
thrust(sol::Solution)
discharge_current(sol::Solution, frame::Integer)
discharge_current(sol::Solution)
anode_eff(sol::Solution, frame::Integer)
anode_eff(sol::Solution)
divergence_eff(sol::Solution, frame::Integer)
divergence_eff(sol::Solution)
ion_current(sol::Solution, frame)
ion_current(sol::Solution)
electron_current(sol::Solution, frame)
electron_current(sol::Solution)
current_eff(sol::Solution, frame)
current_eff(sol::Solution)
mass_eff(sol::Solution, frame)
mass_eff(sol::Solution)
voltage_eff(sol::Solution, frame::Integer)
voltage_eff(sol::Solution) 
```
