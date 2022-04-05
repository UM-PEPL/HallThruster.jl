# Configuration

The `Config` struct contains all of the options you need to run a simulation. On this page, we will explain what options are available and what they do. Note that all arguments must be provided as keywords.

There are five absolutely mandatory arguments. These are:

- `discharge_voltage`: The difference in potential between the anode and cathode, in Volts. This is used to set the left boundary condition. If the cathode potential is zero, then the anode potential is equal to the discharge voltage.
- `thruster`: This is a `Thruster` object containing important geometric and magnetic information about the thruster being simulated. See the page about Thrusters for more.
- `domain`: This is a Tuple containing the locations of the left and right boundaries of the simulation domain, in meters. For instance, if your simulation domain starts at z = 0.0 and is 5 cm long, you would write `domain = (0.0, 0.05)`.
- `anode_mass_flow_rate`: The propellant mass flow rate at the anode, in kg/s
- `initial_condition!`: A function used for initializing the simulation. See the page about [Initialization](initialization.md) for more information.

Aside from these arguments, all others have  default values provided. These are detailed below:

- `ncharge`: Number of charge states to simulate. Defaults to `1`.
- `propellant`: Propellant gas. Defaults to `Xenon`. Other options are described on the Gases and Species page.
- `scheme`: Numerical scheme to employ for integrating the ion equations. This is a `HyperbolicScheme` struct with fields `flux_function`, `limiter`, and `reconstruct`. Defaults to `HyperbolicScheme(flux = rusanov, limiter = minmod, reconstruct = false)`. For more information, see [Fluxes and Numerics](@ref).
- `cathode_potential`: The potential at the right boundary of the simulation. Defaults to `0.0`
- `anode_Te`: The electron temperature at the anode, in eV. Acts as a Dirichlet boundary condition for the energy equation. Defaults to `3.0`.
- `cathode_Te`: The electron temperature at the cathode, in eV. Acts as a Dirichlet boundary condition for the energy equation. Defaults to `3.0`.
- `wall_loss_model`: How radial losses due to sheaths are computed. Defaults to `ConstantSheathPotential(sheath_potential=-20.0, inner_loss_coeff = 1.0, outer_loss_coeff = 1.0)`, which is the loss term from LANDMARK case 1. Other wall loss models are described on the [Wall Loss Models](@ref) page.
- `wall_collision_freq`: Extra "wall collisions" to be added to the total electron momentum transfer collision frequency inside of the channel.  Units of Hz. Defaults to `0.0`.
- `anom_model`: Model for computing the anomalous collision frequency. Defaults to `TwoZoneBohm(1/160, 1/16)`. Further details on the [Anomalous Transport](@ref) page.
- `ionization_model`: Model for ionization. Defaults to `IonizationLookup()`, which uses a lookup table to compute ionization rate coefficients as a function of electron energy. Other options are described on the [Collisions and Reactions](@ref) page.
- `excitation_model`: Model for excitation reactions. Defaults to `ExcitationLookup()`, which uses a lookup table to compute excitation rate coefficients as a function of electron energy.. Other models are described on the [Collisions and Reactions](@ref) page.
- `electron_neutral_model`: Model for elastic scattering collisions between electrons and neutral atoms. Defaults to `ElectronNeutralLookup()`, which uses a lookup table to compute the elastic scattering rate coefficient. Other models are described on the [Collisions and Reactions](@ref) page.
- `electron_ion_collisions`: Whether to include electron-ion collisions. Defaults to `true`. More information on the [Collisions and Reactions](@ref) page.
- `neutral_velocity`: Neutral velocity in m/s. Defaults to `300.0`
- `neutral_temperature`: Neutral temperature in Kelvins. Defaults to `300.0`.
- `ion_temperature`: Ion temperature in Kelvins. Defaults to 100.0
- `electron_pressure_coupled`: Whether to use an electron-pressure-coupled method for the ions.  When using a coupled scheme, ion speed of sound is raised from the thermal speed to the ion acoustic speed, which reduces numerical oscillations around the ion stagnation point. Partial coupling is possible by supplying a value between zero and one, where zero is fully uncoupled and one is fully coupled. Described further in [K. Hara, *Non-oscillatory quasineutral fluid model of cross-field discharge plasmas*, Physics of Plasmas 25, 123508, 2018](https://aip.scitation.org/doi/pdf/10.1063/1.5055750). Defaults to `1.0`.
- `implicit_energy`: The degree to which the energy is solved implicitly. `0.0` is a fully-explicit forward Euler, `0.5` is Crank-Nicholson, and `1.0` is backward Euler. Defaults to `1.0`.
- `min_number_density`: Minimum allowable number density for any species. Defaults to `1e6`
- `min_electron_temperature`: Minimum allowable electron temperature. Defaults to `1.0`.
- `callback`: User-provided callback. This can by any standard callback from `DifferentialEquations.jl`. Defaults to `nothing`.
- `magnetic_field_scale`: Factor by which the magnetic field is increased or decreased compared to the one in the provided `Thruster` struct. Defaults to `1.0`.
- `source_neutrals`: Extra user-provided neutral source term. Can be an arbitrary function, but must take `(U, params, i)` as arguments. Defaults to `Returns(0.0)`. See [User-Provided Source Terms](@ref) for more information.
- `source_ion_continuity`: Vector of extra source terms for ion continuity, one for each charge state. Defaults to `fill(Returns(0.0), ncharge)` . See [User-Provided Source Terms](@ref) for more information.
- `source_ion_momentum`: Vector of extra source terms for ion momentum, one for each charge state. Defaults to `fill(Returns(0.0), ncharge)` . See [User-Provided Source Terms](@ref) for more information.
- `source_potential`: Extra source term for potential equation. Defaults to `Returns(0.0)`. See [User-Provided Source Terms](@ref) for more information.
- `source_electron_energy`: Extra source term for electron energy equation. Defaults to `Returns(0.0)`. See [User-Provided Source Terms](@ref) for more information.
