# Anomalous Transport

HallThruster has a few anomalous transport models built in and allows users to define their own.
Currently, we only support anomalous transport models that are fixed as a function of space or
algebraic models which depend on the local plasma properties, but multi-equation transport
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
\nu_{AN} = \begin{cases}
    \c_1 \omega_{ce} & z < L_{ch} \\
    \c_2 \omega_{ce} & z \ge L_{ch}
\end{case}
```

In the above expression, ``c_1`` and ``c_2`` are tunable coefficients, ``omega_{ce} = e B / m_e`` is the electron cyclotron frequency, and ``L_{ch}`` is the channel length. A `TwoZoneBohm` model is initialized as follows

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
All anomalous transport models in HallThruster.jl are subtypes of `AnomalousTransportModel`.
Two methods must be defined for all such models.

### `initialize_anom!(νan, model, U, params)`
This function takes a pre-allocated vector `νan` of size `length(params.z_cell)`, the model,
the initial simulation state vector `U`, and the simulation `params`, and fills `νan` with
the anomalous collision frequency at the start of the simulation at each cell center.

### `evaluate_anom(U, params, i)`
This function takes the state vector `U`, the simulation `params` and cell number `i` and
returns the anomalous collision frequency for that cell.

## Types of model

### `FixedAnomModel`
These are models which do not change as the simulation progresses. To define a new subtype
of `FixedAnomModel`, you only need to implement `initialize_anom!` as `evaluate_anom` is by
default not called for subtypes of `FixedAnomModel`.

### `ZeroEquationModel`
These depend on the local plasma parameters. By default, `initialize_anom!` for such models
will just call the user-implemented `evaluate_anom` function at each grid location and does
not need to be implemented by the user.

## Custom anomalous transport models
Users of HallThruster may define their own models by defining a custom subtype
of `AnomalousTransportModel`. Suppose we want to implement ``nu_{AN} = \beta \omega_{ce}``
(classic Bohm diffusion). This is a fixed anomalous transport model and does not change
as the simulation progresses. We would first define our type:

```julia
using HallThruster

struct BohmDiffusion <: FixedAnomModel
    β::Float64
end
```

We then need to define the `initialize_anom!` function

```julia
function initialize_anom!(νan, model::BohmDiffusion, U, params)
    for i in 1:length(params.z_cell)
        B = params.cache.B[i]
        ωce = HallThruster.e * B / HallThruster.me
        νan[i] = model.β * ωce
    end
    return νan
end
```

Alternatively, we could define this as a `ZeroEquationModel`. This would be wasteful as
we do not need to recompute the anomalous collision frequency at every iteration, but we
present it for completeness. For `ZeroEquationModel` types, we only need to define
`evaluate_anom` and we do it as so:

```julia
function evaluate_anom(model::BohmDiffusion, U, params, i)
    B = params.cache.B[i]
    ωce = HallThruster.e * B / HallThruster.me
    νan[i] = model.β * ωce
    return νan
end

```

Regardless of which type we choose,, we can use this in a `Config` struct (see [Configuration](@ref))
and the simulation will correctly compute the anomalous transport according to our model.
