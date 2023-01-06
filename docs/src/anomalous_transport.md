# Anomalous Transport

HallThruster has a few anomalous transport models built in and allows users to define their own.
Currently, we only support anomalous transport models that are fixed as a function of space or algebraic models which depend on the local plasma properties, but multi-equation transport
models are planned for the future.

## Built-in Models

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
    \nu_{AN} &= \c_1 \omega_{ce} \quad z < L_{ch} \\
    &= \c_2 \omega_{ce} \quad z \ge L_{ch}
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

We then need to define the model

```julia
function (model::BohmDiffusion)(U, params, i)
    B = params.cache.B[i]
    ωce = HallThruster.e * B / HallThruster.me
    νan = model.β * ωce
    return νan
end
```

We can now set `anom_model = BohmDiffusion` in our config struct (see [Configuration](@ref)) and the simulation will correctly compute the anomalous transport according to our model.
