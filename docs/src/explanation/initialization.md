# Initialization

HallThruster.jl provides sensible defaults for simulation initialization, or allows you to specify your own initial condition.

## Default

The default is `DefaultInitialization()`, which initializes the solution domain as described in the following sections. Below, ``z_0`` and ``z_N`` are `domain[1]` and `domain[2]`, as passed into the `Config` object (see [Configuration](@ref)), ``L_{ch}`` and ``A_{ch}`` are `config.thruster.geometry.channel_length` and `config.thruster.geometry.channel_area`, respectively, and ``\dot{m}`` is `config.anode_mass_flow_rate`.

The `DefaultInitialization` has parameters `max_electron_temperature`, `min_ion_density`, and `max_ion_density`, which can be used to scale the default initialized.
These default to `config.discharge_voltage/10`, `2e17`, and `1e18`, respectively.

### Ion densities

The ion densities are Gaussian with a constant offset and a scaling factor proportional to the mass flow rate and discharge voltage.  For ions with charge 1, the density is
```math
\rho_{i} = m_i \sqrt{\frac{V_d}{300}}\frac{\dot{m}}{5\times10^{-6}}\left(\rho_{min}  + (\rho_{max}-\rho_{min}) \exp\left[-\left(\frac{z - z_0 - L_{ch}/2}{L_{ch}/3}\right)^2\right]\right)
```
For ions with charge `Z`, the density is assumed to scale as
```math
\rho_i|_Z = \frac{\rho_i |_{Z=1}}{Z^2}
```

### Ion velocities

Ions are initialized with the Bohm velocity at the anode. For an ion of charge ``Z``, this is
```math
u_i[1] = -u_{bohm} =- \sqrt{\frac{Z \;e\;T_{eV, anode}}{m_i}}
```


The maximum ion velocity is determined by the discharge voltage ``V_d``:
```math
u_i[\mathrm{end}] = u_{max} = \sqrt{\frac{2 \;Z \;e \; V_d}{m_i}}
```
The initial ion velocity profile between the cathode and the anode is then prescribed as:
```math
u_i(z) = \begin{cases}
	u_{bohm} + \frac{2}{3}(u_{max} - u_{bohm})\left(\frac{z - z_0}{L_{ch}}\right)^2 & z-z_0 < L_{ch} \\
	\frac{1}{3}\left(u_{bohm} + u_{max}\right)\left(1 - \frac{z - z_0 - L_{ch}}{z_N - L_{ch}}\right) + u_{max}\left(\frac{z - z_0 - L_{ch}}{z_N - L_{ch}}\right) & z - z_0 \ge L_{ch}
\end{cases}
```

### Neutral density

The neutral density at the anode is computed in the same way as during a simulation, namely:
```math
\rho_{n, anode} = \frac{\dot{m}}{u_n A_{ch}} - \sum_s \frac{[\rho_{is} u_{is}]_{anode}}{u_n}
```
The density at the cathode is assumed to be 1/100 that at the anode. In the domain, the neutral density has a sigmoid shape:
```math
\rho_n(z) = \frac{1}{2}\left(\rho_{n,anode} + \rho_{n, cathode} + (\rho_{n, cathode} - \rho_{n, anode})\tanh\left(\frac{z - z_0 - L_{ch}/2}{L_{ch} / 24}\right)\right)
```

### Electron energy

The number density is computed from the ion densities. The electron temperature is a Gaussian with height ``V_d / 10`` eV plus a linear baseline to make sure the boundary conditions are satisfied:
```math
T_e(z) = \left(1 - \frac{z - z_0}{z_N - z_0}\right) T_{e, anode} + \left(\frac{z - z_0}{z_N - z_0}\right) T_{e, cathode} + (T_{e,max} - T_{e,min})\exp\left[-\left(\frac{z - z_0 - L_{ch}}{L_{ch}/3}\right)^2\right]
```

### Example

For  a simulation of the SPT-100 with ``V_d``= 500V, three ion charge states, a a mass flow rate of 3 mg/s, an anode electron temperature of 3 eV and a cathode electron temperature of 5 eV, the initial condition looks like:

![](https://github.com/UM-PEPL/HallThruster.jl/blob/main/docs/src/assets/intialization.jpg?raw=true)

## Custom initial conditions

You may define your own initial condition by creating subtypes of `HallThruster.InitialCondition`. Let's say for some reason we wanted to initialize every fluid's density and momentum in every cell to the z-location of its cell center. We might define our initialization as follows:

```jldoctest initialization; output=false
using HallThruster

struct MyInitialCondition <: HallThruster.InitialCondition end;

# output

```

We would then add a method to the `initialize!(U, params, model)` function as follows:

```jldoctest initialization; output=false
import HallThruster.initialize!

function HallThruster.initialize!(params, config, model::MyInitialCondition)
	(;grid, fluid_containers) = params # Pull cell centers locations out of params

    for fluid in fluid_containers.continuity
        fluid.density .= grid.cell_centers
    end

    for fluid in fluid_containers.isothermal
        fluid.density .= grid.cell_centers
        fluid.momentum .= grid.cell_centers
    end

    return nothing
end;

# output

```

We can then test out our initialization function

```jldoctest initialization

config = HallThruster.Config(
    ncharge = 1,
    thruster = HallThruster.SPT_100,
    domain = (0, 0.08),
    anode_mass_flow_rate = 5e-6,
    discharge_voltage = 300.0,
    initial_condition = MyInitialCondition(),
)

simparams = HallThruster.SimParams(
    duration = 1e-3,
    grid = HallThruster.EvenGrid(100)
)

params = HallThruster.setup_simulation(config, simparams)


# the boundary nodes will differ from what we specified due to the enforcement of boundary conditions
params.fluid_containers.continuity[1].density[2:end-1] == params.grid.cell_centers[2:end-1] &&
params.fluid_containers.isothermal[1].density[2:end-1] == params.grid.cell_centers[2:end-1] &&
params.fluid_containers.isothermal[1].momentum[2:end-1] == params.grid.cell_centers[2:end-1]

# output

true
```