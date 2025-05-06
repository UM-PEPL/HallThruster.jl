```@meta
EditURL = "NEWS.md"
```

# Release notes

## v0.18.6

### New features
- Added molecular nitrogen (`HallThruster.MolecularNitrogen`) propellant along with ionization and elastic scattering cross sections.

## v0.18.5
### Bug fixes
- Remove debug print left in prev version

## v0.18.4

### New features
#### Anomalous transport models
A new `ScaledGaussianBohm` transport model has been added.
This is a reparameterized version of the `GaussianBohm` model with non-dimensional parameters of order unity.
See [Anomalous transport](@ref) for more details.


## v0.18.3

### Bug fixes
- Fix indexing error in rate coefficient calculation
- Fix massive divergence angles when ion velocity is zero or negative

## v0.18.2

### New features
#### Anomalous transport models
A new `SimpleLogisticShift` pressure shift model has been added.
This is a simplified version of the `LogisticPressureShift` with one fewer parameter.
See [Anomalous transport](@ref) for more details.

## v0.18.1

### Minor changes
#### Anomalous transport models
The pressure shift now only applies to the coefficients of the transport model as opposed to shifting the entire magnetic field.
This was the pre-0.18.0 behavior.

## v0.18.0

!!! warning "Breaking release"
    v0.18.0 is a breaking release, made as part of our effort to move toward v1.0.0 in the next few months.
    Users may need to update their code to avoid errors.
    Check out the list of changes and removals below to see how to migrate your code.

### Major changes

#### Precompilation and load time
- We have significantly reduced our load times by removing dependencies. We have also added a significant amount of additional precompilation work to reduce the time-to-first-simulation for most users.
- With these improvements, HallThruster.jl may be suitable for use in a scripting workflow, depending on the expected run-time of your simulations.
- For many users and workflows, the total time-to-first-simulation should be just a few seconds.

#### Exports
`Xenon`, `Krypton`, and `time_average` are no longer exported.
Users will need to use the fully-qualified names, i.e. `HallThruster.Xenon`, `HallThruster.Krypton`, and `HallThruster.time_average` to access these variables.
This can be made less verbose by bringing `HallThruster` into scope with a shorter name:
```julia
using HallThurster: HallThruster as het

het.Xenon
```

#### Postprocessing changes
- The `run_simulation` function now takes an optional `postprocess` keyword argument, which takes a `Postprocess` object.
- The `Postprocess` struct allows the user to specify time-averaging, saving, and output info. Any requested postprocessing will be performed automatically after the simulation has finished and output the requested JSON file.
- See [Postprocessing](@ref) for more information

#### JSON IO
JSON input and output has been significantly reworded and enhanced.
See  [Use JSON for input and output](@ref) for more information.
Users of the old JSON interface will need to update their input files to conform to the new format.
Some of the changes include:
- JSON inputs now map directly onto the `Config` structs used by HallThruster.jl elsewhere. The key names are now identical between the JSON input and the `Config` struct.
- JSON outputs now contain all input used to run the simulation, including a `Config` (under key `"config"`) and timestepping info (under key `"simulation"`).
- JSON outputs can be used to restart simulations, replacing the restart capability previously provided by `JLD2`.
- Users can provide postprocessing directives in the JSON input by providing a `"postprocessing"` key in their input, following the [`Postprocessing`](@ref) interface introduced in this update.

### Python interface

We have a new python interface. See [Run a simulation from python](@ref) for more.

### Interface changes

#### Configuration changes
- The following keys have been renamed to better express the expected units.
    - `anode_Te` -> `anode_Tev`
    - `cathode_Te` -> `cathode_Tev`
    - `ion_temperature` -> `ion_temperature_K`
    - `neutral_temperature` -> `neutral_temperature_K`
    - `background_neutral_temperature` -> `background_temperature_K`
    - `background_pressure` -> `background_pressure_Torr`
    - `cathode_potential` -> `cathode_coupling_voltage`
- The following keys have been removed.
    - `min_electron_temperature` (now set to the `min(anode_Tev, cathode_Tev)`)
    - `min_number_density` (now set to `1e6`)
- The following keys have had their options changed
    - `ionization_model`: Formerly took a value of type `IonizationModel`, now takes a `Symbol` (one of `:Lookup` or `:Landmark`)
    - `excitation_model`: Formerly took a value of type `ExcitationModel`, now takes a `Symbol` (one of `:Lookup` or `:Landmark`)
    - `electron_neutral_model`: Formerly took a value of `ElectronNeutralModel`, now takes a `Symbol` (one of `None`, `:Lookup` or `:Landmark`)
- See [Configuration](@ref) for a full list of options and default values.

#### Simulation changes
- A new struct `SimParams` has been created to hold timestepping and grid generation options
- `run_simulation` can now be called as `run_simulation(config::Config, simparams::SimParams)`.
- See [Simulations](reference/simparams.md) for more information.

#### Solution changes
- `u` has been removed from `Solution`
- `Solution.savevals` has been renamed to `frames`
- `Solution.config` contains the `Config` used to run the simulation
- `Solution.error` captures any errors that are produced during a run
- Retcode options have changed. The options are now `:success`, `:failure`, and `:error`. `:success` denotes a successful simulation, as before. `:failure` indicates a numerical instability, and `:error` captures any cases in which the simulation would have thrown an exception. The error text is captured in `Solution.error`.

#### Wall loss model changes
- The parameter `α` of the `WallSheath` model has been renamed to `loss_scale` and has a default value of 1.0.
- See [Wall Loss Models](@ref) for more information.

#### Grid generation changes
- `UnevenGrid` no longer takes a custom density function as an input.
- See [Grid generation](explanation/grid_generation.md) for more information.

#### Anomalous transport model changes
- `ShiftedTwoZoneBohm`, `ShiftedGaussianBohm`, and `ShiftedMultiBohm` have been removed.
- The model `LogisticPressureShift` has been added. This takes another anom model as a parameter and applies the same pressure-dependent shift as seen in these models.
- To migrate to the new system, replace `ShiftedGaussianBohm` or similar with `LogisticPressureShift(GaussianBohm(...), ...)`.
- The `dz` and `z0` parameters of `LogisticPressureShift` are now in terms of channel length rather than meters
- Added the `GaussianBohm` model
- Renamed the parameters of the `ShiftedGaussianBohm`/`GaussianBohm` model. The new parameters are:
    - `hall_min`: the minimum Hall parameter
    - `hall_max`: the maximum Hall parameter
    - `width`: the width of the Gaussian subtraction from the `Bohm`-like transport model.
    - `center`: the center of the Gaussian subtraction
- See [Anomalous transport](@ref) for a full listing of anomalous transport options.

#### Thruster changes
- `Thruster.magnetic_field` now takes a `MagneticField` object, which contains a filename, z-coordinates, and magnetic field values. Arbitrary functions are no longer supported. If only a filename is provided, the constructor will load the magnetic field from the provided file.
- `B_field_SPT_100` has been renamed to `spt100_analytic_field`
- `geometry_SPT_100` has been renamed to `spt100_geometry`
- See [Thrusters](@ref) for more information.

### Removals
- The `IonizationModel`, `ExcitationModel`, and `ElectronNeutralModel` interfaces have been removed in favor of a simpler `Symbol`-based lookup system.
- The `GKElectronNeutral` model has been removed entirely.
- Serialization using JLD2 has been removed in favor of the JSON-based system described above.
- Plotting capabilities are no longer provided by HallThruster.jl directly.
