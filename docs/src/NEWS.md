# HallThuster.jl v0.18.0

This is a breaking release, as part of our effort to move toward v1.0.0.
Users may need to update their code to avoid errors. 
Check out the list of changes and removals below to see how to migrate your code.

## Major changes

### Precompilation and load time
- We have significantly reduced our load times by removing dependencies. We have also added a significant amount of additional precompilation work to reduce the time-to-first-simulation for most users.
- With these improvements, HallThruster.jl may be suitable for use in a scripting workflow, depending on the expected run-time of your simulations.
- For many users and workflows, the total time-to-first-simulation should be just a few seconds.

### Postprocessing
- The `run_simulation` function now takes an optional `postprocess` keyword argument, which takes a `Postprocess` object.
- The `Postprocess` struct allows the user to specify time-averaging, saving, and output info. Any requested postprocessing will be performed automatically after the simulation has finished and output the requested JSON file.

### JSON
JSON input and output has been significantly reworded and enhanced. 
See the page on [JSON input and output](json.md) for more information.
Users of the JSON interface will need to update their input files to conform to the new format.
Some of the changes include:
- JSON inputs now map directly onto the `Config` structs used by HallThruster.jl elsewhere. The key names are now identical between the JSON input and the `Config` struct.
- JSON outputs now contain all input used to run the simulation, including a `Config` (under key `"config"`) and timestepping info (under key `"simulation"`).
- JSON outputs can be used to restart simulations, replacing the restart capability previously provided by `JLD2`.
- Users can provide postprocessing directives in the JSON input by providing a `"postprocessing"` key in their input, following the [`Postprocessing`](postprocessing.md) interface introduced in this update.

## Interface changes

### Configuration
- The following keys have been renamed to better express the expected units. See [Configuration](configuration.md) for a full list of options and default values.
    - `anode_Te` -> `anode_Tev`
    - `cathode_Te` -> `cathode_Tev`
    - `ion_temperature` -> `ion_temperature_K`
    - `neutral_temperature` -> `neutral_temperature_K`
    - `background_neutral_temperature` -> `background_temperature_K`
    - `background_pressure` -> `background_pressure_Torr`
- The following keys have been removed.
    - `min_electron_temperature` (now set to the `min(anode_Tev, cathode_Tev)`)
    - `min_number_density` (now set to `1e6`)
- The following keys have had their options changed
    - `ionization_model`: Formerly took a value of type `IonizationModel`, now takes a `Symbol` (one of `:Lookup` or `:Landmark`)
    - `excitation_model`: Formerly took a value of type `ExcitationModel`, now takes a `Symbol` (one of `:Lookup` or `:Landmark`)
    - `electron_neutral_model`: Formerly took a value of `ElectronNeutralModel`, now takes a `Symbol` (one of `None`, `:Lookup` or `:Landmark`)

### Wall loss models
- The parameter `α` of the `WallSheath` model has been renamed to `loss_scale` and has a default value of 1.0.

### Grid generation
- `UnevenGrid` no longer takes a custom density function as an input.

### Anomalous transport models
- `ShiftedTwoZoneBohm`, `ShiftedGaussianBohm`, and `ShiftedMultiBohm` have been removed.
- The model `LogisticPressureShift` has been added. This takes another anom model as a parameter and applies the same pressure-dependent shift as seen in these models.
- To migrate to the new system, replace `ShiftedGaussianBohm` or similar with `LogisticPressureShift(GaussianBohm(...), ...)`.
- Added the `GaussianBohm` model
- Renamed the parameters of the `ShiftedGaussianBohm`/`GaussianBohm` model. The new parameters are:
    - `hall_min`: the minimum Hall parameter
    - `hall_max`: the maximum Hall parameter
    - `width`: the width of the Gaussian subtraction from the `Bohm`-like transport model.
    - `center`: the center of the Gaussian subtraction
- See [Anomalous Transport][anomalous_transport.md] for a full listing of anomalous transport options.

### Thrusters
- `Thruster.magnetic_field` now takes a `MagneticField` object, which contains a filename, z-coordinates, and magnetic field values. Arbitrary functions are no longer supported. If only a filename is provided, the constructor will load the magnetic field from the provided file.
- `B_field_SPT_100` has been renamed to `spt100_analytic_field`
- `geometry_SPT_100` has been renamed to `spt100_geometry`

## Removals
- The `IonizationModel`, `ExcitationModel`, and `ElectronNeutralModel` interfaces have been removed in favor of a simpler `Symbol`-based lookup system.
- The `GKElectronNeutral` model has been removed entirely.
- Serialization using JLD2 has been removed in favor of the JSON-based system described above.
- Plotting capabilities are no longer provided by HallThruster.jl directly. 

