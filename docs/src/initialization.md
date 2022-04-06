# Initialization

Hall2De provides sensible defaults for simulation initialization, or allows you to specify your own initial condition.

## Default

The default is `DefaultInitialization()`, which initializes the solution domain as described in the following sections. Below, $z_0$ and $z_N$ are `domain[1]` and `domain[2]`, as passed into the `Config` object (see [Configuration](@ref)), $L_{ch}$ and $A_{ch}$ are `config.thruster.geometry.channel_length` and `config.thruster.geometry.channel_area`, respectively, and $\dot{m}$ is `config.anode_mass_flow_rate`.

### Plasma density

### Ion velocities

Ions are initialized with the Bohm velocity at the anode. For an ion of charge ``Z``, this is
$$
u_i[1] = -u_{bohm} =- \sqrt{\frac{Z \;e\;T_{eV, anode}}{m_i}}
$$


The maximum ion velocity is determined by the discharge voltage ``V_d``:
$$
u_i[\mathrm{end}] = u_{max} = \sqrt{\frac{2 \;Z \;e \; V_d}{m_i}}
$$
The initial ion velocity profile between the cathode and the anode is then prescribed as:
$$
u_i(z) = \begin{cases}
	u_{bohm} + \frac{2}{3}(u_{max} - u_{bohm})\left(\frac{z - z_0}{L_{ch}}\right)^2 & z-z_0 < L_{ch} \\
	\frac{1}{3}\left(u_{bohm} + u_{max}\right)\left(1 - \frac{z - z_0 - L_{ch}}{z_N - L_{ch}}\right) + u_{max}\left(\frac{z - z_0 - L_{ch}}{z_N - L_{ch}}\right) & z - z_0 \ge L_{ch}
\end{cases}
$$

### Neutral density

The neutral density at the anode is computed in the same way as during a simulation, namely: 
$$
\rho_{n, anode} = \frac{\dot{m}}{u_n A_{ch}} - \sum_s \frac{[\rho_i u_i]_{anode}}{u_n}
$$
The density at the cathode is assumed to be 1/100 that at the anode. In the domain, the neutral density has a sigmoid shape:
$$
\rho_n(z) = 
$$




### Electron energy



### Example

For  a simulation of the SPT-100 with $V_d$= 500V, three ion charge states, a a mass flow rate of 3 mg/s, an anode electron temperature of 3 V and a cathode electron temperature of 5 V, the initial condition looks like:

![](./assets/initialization.svg)

## Custom initial conditions

You may define your own initial condition by creating subtypes of `HallThruster.InitialCondition`. Let's say for some reason we wanted to initialize every state variable in every cell to the z-location of its cell center. We might define our initialization as follows:

```jldoctest initialization; output=false
using HallThruster

struct MyInitialCondition <: HallThruster.InitialCondition end;

# output

```

We would then add a method to the `initialize!(U, params, model)` function as follows:

```jldoctest initialization; output=false
import HallThruster.initialize!

function HallThruster.initialize!(U, params, model::MyInitialCondition)
	(;z_cell) = params # Pull cell centers locations out of params
    nvars = size(U, 1)
    for (i, z) in enumerate(z_cell)
       	for j in 1:nvars
           	U[j, i] = z_cell[i]
        end
    end
    return U # optional. Since U is modified. the return value is never used, but by Julia convention we also return the mutated object.
end;

# output

```

We can check the behavior of our new function:

```jldoctest initialization
# Dummy config and params
ncells = 100
nvars = 4
config = (;initial_condition = MyInitialCondition())
z_cell = LinRange(0, 0.05, ncells)
U = zeros(nvars, ncells)
params = (;config, z_cell)

# Method of initialize! which dispatches to initialize!(U, params, config.initial_condition)
# This is what HallThruster.jl calls when initializing a simulation
HallThruster.initialize!(U, params)

U[1, :] == U[2, :] == U[3, :] == U[4, :] == collect(z_cell)

# output

true
```

