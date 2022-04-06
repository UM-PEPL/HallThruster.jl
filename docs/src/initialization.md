# Initialization

Hall2De provides sensible defaults for simulation initialization, or allows you to specify your own initial condition.

## Default

The default is `DefaultInitialization()`, which initializes the solution domain in the following way

### Neutral density



### Plasma density



### Ion velocities



### Electron energy



### Example



For  a simulation of the SPT-100 with $V_d$= 500V, three ion charge states, a a mass flow rate of 3 mg/s, an anode electron temperature of 3 V and a cathode electron temperature of 5 V, the initial condition looks like:

![initialization_img](./assets/initialization.svg)

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

