# Anomalous Transport

HallThruster has a few anomalous transport models built in and allows users to define their own. This page describes these models and the process by which algebraic and multi-equation transport
models can be added by the user.

!!! warning "Interface not finalized"
    The AnomalousTransportModel interface is not yet finalized and subject to revision.
    Keep this in mind when using this feature.

## Built-in Models

!!! warning Incomplete documentation
Some models (those relating to pressure-dependent effects) have not been finalized and are not represented in this documentation.
Please see the [source code](https://github.com/UM-PEPL/HallThruster.jl/blob/main/src/collisions/anomalous.jl) for a complete listing.

### `NoAnom()`

Model for no anomalous transport (anomalous collision frequency  = 0).

```julia
anom_model = NoAnom()
```

### `Bohm(c)`

Model where the anomalous collision frequency scales with the electron cyclotron frequency ωce times some scaling factor c

### `TwoZoneBohm(c1, c2)`

HallThruster's default anomalous transport option. This is a standard model of anomalous transport frequently used in Hall thruster simulations. The anomalous collision frequency is defined as

```math
\begin{aligned}
    \nu_{AN} &= c_1 \omega_{ce} \quad z < L_{ch} \\
    &= c_2 \omega_{ce} \quad z > L_{ch}
\end{aligned}
```

In the above expression, ``c_1`` and ``c_2`` are tunable coefficients, ``\omega_{ce} = e B / m_e`` is the electron cyclotron frequency, and ``L_{ch}`` is the channel length. A `TwoZoneBohm` model is initialized as follows

```julia
anom_model = TwoZoneBohm(c1, c2)
```

The transition between the zones is determined by the user-provided transition function.
This defaults to a step function.

### `MultiLogBohm(z, c)`

Model similar to that employed in Hall2De, where the mobility is Bohm-like (i.e. `νan(z) = c(z) * ωce(z)`) and z is in meters.

The function `c(z)` is defined by a sequence of nodes `(z, c)` provided by the user. At `z = z[1]`, `c(z) = c[1]`, and so forth.

At `z[i] < z < z[i+1]`, `log(c)` is defined by linearly interpolating between `log(c[i])` and `log(c[i+1])`.

For `z < z[1]`, `c = c[1]` and for `z > z[end]`, `c(z) = c[end]`.

The user may also provide a single array of `[z[1], z[2], ..., z[end], c[1], c[2], ..., c[end]]`. The number of `z `values must be equal to the number of c values.

## The `AnomalousTransportModel` interface
Currently, HallThruster.jl expects all models to be written to be callable structs, taking arguments `U, params, i`, where `U` is the system state vector, `params` are the simulation parameters (including the cache of all variables), and `i` is the index of the cell.

## Custom anomalous transport models
Users of HallThruster may define their own models by defining a custom subtype
of `AnomalousTransportModel`. Suppose we want to implement ``nu_{AN} = \beta \omega_{ce}``
(classic Bohm diffusion). This is a fixed anomalous transport model and does not change
as the simulation progresses. We would first define our type:

```julia
using HallThruster

struct BohmDiffusion <: AnomalousTransportModel
    β::Float64
end
```

We then need to define the function which computes the anomalous transport in each cell.
This is a mutating function which takes two arguments: the vector of anomalous collision
frequency values to be updated, and the solver params.

```julia
function (model::BohmDiffusion)(νan, params)

    e = HallThruster.e
    me = HallThruster.me
    B = params.cache.B

    for i in eachindex(νan)
        ωce = e * B[i] / me
        νan[i] = model.β * ωce
    end

    return νan
end
```

We can now set `anom_model = BohmDiffusion` in our config struct (see [Configuration](@ref)) and the simulation will correctly compute the anomalous transport according to our model.

## More complex anomalous transport models

Up until now, we have only defined algebraic models of anomalous electron transport.
For more high-fidelity models, we might need to solve multiple partial differential equations.
Even if we don't want to do that, we still might want to compute other anomalous transport-
related quantities at the same time that we update the anomalous transport. Lets see how
we might do this.

As a first example, let's compute an anomalous transport that depends on energy density of
electrostatic waves in the plasma. This model derives from Lafleur, Chabert, and Balruud (2016)
and has the following form:

```math
\begin{aligned}
    \nu_{AN} = K \frac{\nabla \cdot \mathbf{u}_i W}{m_e n_e c_s v_{de}}
\end{aligned}
```

In this expression, ``K`` is a tunable coefficient, ``u_i`` is the ion velocity,
``W = n_e k_B T_e`` is the wave energy density, ``m_e`` is the electron mass, ``n_e`` is the electron number density, ``c_s`` is the ion sound speed, and ``v_{de}`` is the electron azimuthal drift speed. Let's say we want to save the wave energy density in addition to
the anomalous collision frequency. We begin by defining the model:

```julia
    struct LafleurModel <: AnomalousTransportModel
        K::Float64
    end
```

Next, we add a method to the `num_anom_variables` function. Since we want to save the
wave energy density, we need 1 additional anomalous transport variable.

```julia
    num_anom_variables(::LafleurModel) = 1
```

Now, we define the behavior of the model in a function. Since the model is based
on assuming the wave energy convects with the ions, we will use upwind differencing for the gradient.

```julia
function (model::LafleurModel)(νan, params)

    (;config, cache) = params
    mi = config.propellant.m
    K = model.K
    e = HallThruster.e
    me = HallThruster.me
    (;ne, Tev, ue, ui, νe, anom_variables) = cache

    ncells = length(νan)

    for i in 2:ncells-1

        W = e * ne[i] * Tev[i]
        Hall_param = e * B[i] / me / νe[i]
        vde = Hall_param * ue[i]
        cs = sqrt(e * Tev[i] / mi)

        # Upwind differencing of gradient term
        if ui > 0
            dz = params.z_cell[i] - params.z_cell[i-1]
            W_left = e * cache.ne[i-1] * cache.Te[i-1]
            grad_ui_W = (ui[1, i] * W[i] - cache.ui[1, i-1] * W_left) / dz
        else
            dz = params.z_cell[i+1] - params.z_cell[i]
            W_right = e * cache.ne[i+1] * cache.Te[i+1]
            grad_ui_W = (cache.ui[1, i+1] * W_right - ui[i] * W[i]) / dz
        end

        # Save W to cache.anom_variables[1]
        anom_variables[1][i] = W

        # Return anomalous collision frequency
        νan[i] = abs(K * grad_ui_W / (me * cs * vde * ne))
    end

    # Neumann BC anomalous transport
    anom_variables[1][1] = anom_variables[1][2]
    anom_variables[1][end] = anom_variables[1][end-1]
    νan[1] = νan[2]
    νan[end] = νan[end-1]

    return νan
end
```

The saved value of the wave energy density can then be recovered as

```julia
    solution.savevals[frame].anom_variables[1]
```

where `solution` is the `Solution` object resulting from a call to `run_simulation`

## Solving arbitrary PDEs using the AnomalousTransportModel interface

We can use this basic interface to solve PDEs within HallThruster.jl. We demonstrate this by solving the scalar advection equation using first-order upwind differencing in space and forward
Euler integration in time, with periodic boundary conditions. The scalar advection equation is given by:

```math
    \begin{aligned}
        \frac{\partial u}{\partial t} + a \frac{\partial u}{\partial x} = 0
    \end{aligned}
```

To begin, we define the model struct and the number of variables we need.

```julia
using HallThruster, Plots

struct ScalarAdvection{F} <: HallThruster.AnomalousTransportModel
    advection_velocity::Float64  # Advection advection_velocity
    initializer::F  # Initialization function
end

# Save two auxilliary variables
#   1) Advected quantity u
#   2) gradient of advected quantity (du/dz)
HallThruster.num_anom_variables(::ScalarAdvection) = 2
```

Next, we define the model function:

```julia
function (model::ScalarAdvection)(νan, params)

    ncells = length(νan)

    # Extract variables from params
    cache = params.cache
    z = params.z_cell
    u = cache.anom_variables[1]
    du_dz = cache.anom_variables[2]
    a = model.advection_velocity
    dt = params.dt

    if params.iteration[] < 1
        # Initialize
        model.initializer(u, z)
    else
        for i in eachindex(νan)
            # Setup for periodic boundary conditions
            if i == 1
                i_minus_1 = ncells
                dz_minus = z[2] - z[1]
            else
                i_minus_1 = i - 1
                dz_minus = z[i] - z[i-1]
            end

            if i == ncells
                i_plus_1 = 1
                dz_plus = z[2] - z[1]
            else
                i_plus_1 = i + 1
                dz_plus = z[i+1] - z[i]
            end

             # Update gradient
            if a > 0
                du_dz[i] = (u[i] - u[i_minus_1]) / dz_minus
            else
                du_dz[i] = (u[i_plus_1] - u[i]) / dz_plus
            end
        end

        # Update advected quantity
        for i in eachindex(νan)
            u[i] -= a * du_dz[i] * dt
        end
    end

    # Return a two-zone bohm anomalous transport result,
    # since we don't really care about the anomalous transport.
    return HallThruster.TwoZoneBohm(1/160, 1/16)(νan, params)
end
```

And that's it! Now all there is to do is define our simulation parameters and run!

```julia

# Set up config
advection_velocity = 1e5
L = 0.08

config = HallThruster.Config(
    domain = (0.0, L),
    anom_model = ScalarAdvection(advection_velocity, initializer),
    thruster = HallThruster.SPT_100,
    discharge_voltage = 300.0,
    anode_mass_flow_rate = 5e-6
)

# Define dt such that CFL condition is obeyed
ncells = 200
dx = L / ncells
CFL = 0.9
dt = min(1e-8, dx * CFL / advection_velocity)
nsteps = 1000

# Run simulation
solution = HallThruster.run_simulation(
    config;
    ncells,
    dt,
    duration = nsteps * dt,
    nsave = nsteps
)

# Extract variables from solution
z = solution.params.z_cell
u = [saveval.anom_variables[1] for saveval in solution.savevals]
```

We can now visualize the results to make sure everything worked well.

```julia
using Plots

# Time needed to transit the domain
t_transit = L / advection_velocity

# Number of periods
num_periods = floor(Int, nsteps * dt / t_transit)

# Plot results
p = plot(; framestyle = :box, xlabel = "x", ylabel = "u", title = "First order upwind for scalar advection")
for i in 0:num_periods
    index = round(Int, i * t_transit / dt) + 1
    plot!(
        p, z, u[index], label = "After $i periods",
        linecolor = cgrad(:turbo, num_periods + 1, categorical = true)[i+1]
    )
end
display(p)
```

![](https://github.com/UM-PEPL/HallThruster.jl/blob/main/docs/src/assets/scalar_advection.png?raw=true)

This looks correct! In this case, we haven't coupled our PDE solution to the anomalous transport, but one could easily do this. In the same way, systems of two, three, or more coupled PDEs can be solved and related to the anomalous collision frequency.

The full script is reproduced below:

```julia
using HallThruster

struct ScalarAdvection{F} <: HallThruster.AnomalousTransportModel
    advection_velocity::Float64  # Advection advection_velocity
    initializer::F  # Initialization function
end

# Save two auxilliary variables
#   1) Advected quantity u
#   2) gradient of advected quantity (du/dz)
HallThruster.num_anom_variables(::ScalarAdvection) = 2

function (model::ScalarAdvection)(νan, params)

    ncells = length(νan)

    # Extract variables from params
    cache = params.cache
    z = params.z_cell
    u = cache.anom_variables[1]
    du_dz = cache.anom_variables[2]
    a = model.advection_velocity
    dt = params.dt

    if params.iteration[] < 1
        # Initialize
        model.initializer(u, z)
    else
        for i in eachindex(νan)
            # Setup for periodic boundary conditions
            if i == 1
                i_minus_1 = ncells
                dz_minus = z[2] - z[1]
            else
                i_minus_1 = i - 1
                dz_minus = z[i] - z[i-1]
            end

            if i == ncells
                i_plus_1 = 1
                dz_plus = z[2] - z[1]
            else
                i_plus_1 = i + 1
                dz_plus = z[i+1] - z[i]
            end

             # Update gradient
            if a > 0
                du_dz[i] = (u[i] - u[i_minus_1]) / dz_minus
            else
                du_dz[i] = (u[i_plus_1] - u[i]) / dz_plus
            end
        end

        # Update advected quantity
        for i in eachindex(νan)
            u[i] -= a * du_dz[i] * dt
        end
    end

    # Return a two-zone bohm anomalous transport result,
    # since we don't really care about the anomalous transport.
    return HallThruster.TwoZoneBohm(1/160, 1/16)(νan, params)
end


# Define initializer function, which is a step function
# between z = 0.01 and z = 0.02
function initializer(u, z)
    for i in eachindex(u)
        if 0.01 < z[i] < 0.02
            u[i] = 1.0
        else
            u[i] = 0.0
        end
    end
    return u
end

# Set up config
advection_velocity = 1e5
L = 0.08

config = HallThruster.Config(
    domain = (0.0, L),
    anom_model = ScalarAdvection(advection_velocity, initializer),
    thruster = HallThruster.SPT_100,
    discharge_voltage = 300.0,
    anode_mass_flow_rate = 5e-6
)

# Define dt such that CFL condition is obeyed
ncells = 200
dx = L / ncells
CFL = 0.9
dt = min(1e-8, dx * CFL / advection_velocity)
nsteps = 1000

# Run simulation
solution = HallThruster.run_simulation(
    config;
    ncells,
    dt,
    duration = nsteps * dt,
    nsave = nsteps
)

# Extract variables from solution
z = solution.params.z_cell
u = [saveval.anom_variables[1] for saveval in solution.savevals]

using Plots

# Time needed to transit the domain
t_transit = L / advection_velocity

# Number of periods
num_periods = floor(Int, nsteps * dt / t_transit)

# Plot results
p = plot(; framestyle = :box, xlabel = "x", ylabel = "u", title = "First order upwind for scalar advection")
for i in 0:num_periods
    index = round(Int, i * t_transit / dt) + 1
    plot!(
        p, z, u[index], label = "After $i periods",
        linecolor = cgrad(:turbo, num_periods + 1, categorical = true)[i+1]
    )
end
display(p)
```
