# Anomalous Transport

HallThruster has two anomalous transport models built in and allows users to define their own. Currently, we only support algebraic (zero-equation) models, but support for models involving additional PDEs will be added in the future.

## Built-in Models

### `TwoZoneBohm`

HallThruster's default anomalous transport option. This is a standard model of anomalous transport frequently used in Hall thruster simulations. The anomalous collision frequency is defined as

```math
\nu_{AN} = \begin{cases}
    \beta_1 \omega_{ce} & z < L_{ch} \\
    \beta_2 \omega_{ce} & z \ge L_{ch}
\end{case}
```

In the above expression, ``\beta_1`` and ``beta_2`` are tunable coefficients, ``omega_{ce} = e B / m_e`` is the electron cyclotron frequency, and ``L_{ch}`` is the channel length. A `TwoZoneBohm` model is initialized as follows

```julia
anom_model = TwoZoneBohm(β1, β2)
```

### `NoAnom`

Model for no anomalous transport (anomalous collision frequency  = 0).

## Custom anomalous transport models

Users of HallThruster may define their own models by defining a custom subtype of `ZeroEquationModel`. Let's say we want to implement ``nu_{AN} = \beta \omega_{ce}`` (classic Bohm diffusion). We would first define our type:

```julia
using HallThruster

struct BohmDiffusion <: HallThruster.ZeroEquationModel
    β::Float64
end
```

We then need to make our struct callable. It must take a signature of `(U, params, i)`, where `U` is the state matrix, `params` is a NamedTuple containing all simulation parameters and a cache of variables, and `i` is the cell index.

```julia
function (b::BohmDiffusion)(U, params, i)
    B = params.cache.B[i]
    ωce = HallThruster.e * B / HallThruster.me
    νan = β * ωce
    return νan
end
```

Now, we can use this in a `Config` struct (see [Configuration](@ref)) and the simulation will correctly compute the anomalous transport according to our model.
