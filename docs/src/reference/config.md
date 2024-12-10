# Configuration

The `Config` struct contains all of the options you need specify the physics, geometry, and numerics of a simulation.
On this page, we will explain what options are available and what they do.
As in all parts of the code, dimensional quantities are SI unless explicitly noted, but units may be provided using [`Unitful`](https://github.com/PainterQubits/Unitful.jl) or [`DynamicQuantities`](https://github.com/SymbolicML/DynamicQuantities.jl)

## Mandatory arguments

There are four mandatory arguments.

- `discharge_voltage` The difference in potential between the anode and cathode, in Volts. This is used to set the left boundary condition for the electrostatic potential.
- `thruster` This is a `Thruster` object containing important geometric and magnetic information about the thruster being simulated. See [Thrusters](thrusters.md) for more information.
- `domain`  This is a `Tuple` containing the locations of the left and right boundaries of the simulation domain, in meters. For instance, if your simulation domain starts at z = 0.0 and is 5 cm long, you would write `domain = (0.0, 0.05)`.
- `anode_mass_flow_rate` The propellant mass flow rate at the anode, in kg/s

## Optional arguments

These arguments do not need to be provided, and have sensible default values provided.

- `ncharge`  Maximum ion charge state
    - default: 1.
- `propellant`  A `Gas`. See [Propellants](propellants.md) for more
    - default: `Xenon`.
- `cathode_coupling_voltage`  The potential at the right boundary of the simulation
    - default: `0.0`
- `anode_boundary_condition` Can be either `:sheath` or `:dirichlet`.
    - default: `:sheath`.
- `cathode_Tev`  The electron temperature at the cathode, in eV.
    - default: `2.0`.
- `anode_Tev`  The electron temperature at the anode, in eV if `anode_boundary_condition == :dirichlet`.
    - default: `2.0`.
- `anom_model`  Model for computing the anomalous collision frequency. See [Anomalous Transport](../reference/anomalous_transport.md) for more info.
    - default: `TwoZoneBohm(1/160, 1/16)`.
- `wall_loss_model`  How radial losses due to sheaths are computed. Other wall models are described on the [Wall Loss Models](wall_loss_models.md) page.
    - default: `WallSheath(BNSiO2, 1.0)`.
- `conductivity_model`  Model for the cross-fieldelectron thermal conductivity. See[Electron Thermal Conductivity](../reference/electron_thermal_conductivity.md) for more.
    - default: `Mitchner()`
- `electron_ion_collisions`  Whether to include electron-ion collisions. See [Collisions and Reactions](@ref) for more.
    - default: `true`.
- `neutral_velocity`  Neutral velocity in m/s.
    - default: `300.0`, or if `neutral_temperature` is set, that parameter is used to compute the velocity using a one-sided maxwellian flux approximation.
- `neutral_temperature_K`  Neutral temperature in Kelvins.
    - default: `500.0`.
- `ion_temperature_K`  Ion temperature in Kelvins.
    - default: 1000.0
- `ion_wall_losses`  Whether we model ion losses to the walls.
    - default: `false`.
- `background_pressure`  The pressure of the background neutrals, in Pascals. These background neutrals are injected at the anode to simulate the ingestion of facility neutrals.
    - default: `0.0`
- `background_neutral_temperature`  The temperature of the background neutrals, in K.
    - default: `150.0`.
- `solve_plume`  Whether quasi-1D beam expansion should be modelled outside of the channel. See [Quasi-1D plume model](../explanation/plume.md) for more.
    - default: `false`
- `apply_thrust_divergence_correction`  Whether the thrust output by HallThruster.jl should include a divergence correction factor of `cos(divergence_angle)`.
    - default: `false`.
- `electron_plume_loss_scale`  The degree to which radial electron losses are applied in the plume. See [Wall Loss Models](@ref) for more information.
    - default: 1.
- `scheme`  Numerical scheme to employ for integrating the ion equations. This is a `HyperbolicScheme` struct with fields `flux_function`, `limiter`, and `reconstruct`. See [Schemes](../reference/schemes.md) for more info.
    - default: `HyperbolicScheme(flux_function = rusanov, limiter = van_leer, reconstruct = true)`.
- `transition_length`: Distance over which the transition between inside and outside the channel is smoothed. Affects wall losses as well as two-zone Bohm-like transport models.
    - default: `0.1 * thruster.geometry.channel_length`
- `implicit_energy`  The degree to which the energy is solved implicitly. `0.0` is a fully-explicit forward Euler, `0.5` is Crank-Nicholson, and `1.0` is backward Euler.
    - default: `1.0`.
- `magnetic_field_scale`  Factor by which the magnetic field is increased or decreased compared to the one in the provided `Thruster` struct.
    - default: `1.0`.
- `initial_condition`  An `InitialCondition`; see [Initialization](../explanation/initialization.md) for more information.
    - default: `DefaultInitialization()`.
- `anom_smoothing_iters`  How many times to smooth the anomalous transport profile. Only useful for transport models that depend on the plasma properties.
    - default: `0`

## Verification and validation options
These options are used in code benchmarking and verification and are not usually needed by end-users.
See [Verification and validation](../explanation/verification.md) for an explanation of how we use these to verify the accuracy of the code.

- `ionization_model` Model for ionization reactions.
    - default: `:Lookup`.
- `excitation_model` Model for excitation reactions.
    - default: `:Lookup`. 
- `electron_neutral_model`  Model for elastic scattering collisions between electrons and neutral atoms.
    - default: `:Lookup`.
- `LANDMARK`  Whether we are using the physics model from the LANDMARK benchmark. This affects whether certain terms are included in the equations, such as electron and heavy species momentum transfer due to ionization and the form of the electron thermal conductivity.
    - default: `false`.
- `source_neutrals`  Extra user-provided neutral source term.
    - default: `nothing` 
- `source_ion_continuity`  Vector of extra source terms for ion continuity, one for each charge state.
    - default: `nothing`.
- `source_ion_momentum`  Vector of extra source terms for ion momentum, one for each charge state.
    - default: `nothing`.
- `source_electron_energy`  Extra source term for electron energy equation.
    - default: `nothing`. 
