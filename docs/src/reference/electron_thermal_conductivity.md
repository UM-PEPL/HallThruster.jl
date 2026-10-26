# Electron Thermal Conductivity

HallThruster has a few electron thermal conductivity models built in and allows users to define their own. This page describes these models and the process by which algebraic transport
models can be added by the user.

## Built-in Models

### `Braginskii`

This form of the thermal conductivity follows the result of S. I. Braginskii, in Reviews of Plasma Physics, edited by M. A. Leontovich (Consultants Bureau, New York, 1965), Vol. 1, p. 205.:

```math
\begin{aligned}
    \kappa_{e\perp} & \approx C \frac{neTe (eV) \nu}{m_e \omega_{ce}^2} \\
    \nu = \nu_{ei} + \nu_{en} + \alpha \nu_{AN} \\
\end{aligned}

Where C is a constant that is based on the value of the effective charge for multiple charge states and Table 1 of the Braginskii reference. 
The coefficient $\alpha$ ranges from zero to one and scales the amount that anomalous collisions contribute to heat flux.
This defaults to 1.0.
A `Braginskii` model can be initialized in the user-defined configuration as:

```julia
conductivity_model = Braginskii(anom_scale = 1.0)
```

### `Mitchner()`

This model employs the form described in Mitchner & Kruger "Partially Ionized Gases" (1973), given below

```math
\begin{aligned}
    \kappa_{e\perp} & \approx \frac{1}{1+\Omega_H^2} \frac{2.4}{1+ \frac{\nu_{ie}}{\sqrt{2}\nu}} \frac{neTe (eV)}{\nu m_e} \\
    \nu = \nu_{ei} + \nu_{en} + \alpha \nu_{AN} \\
\end{aligned}
```

As in the Braginskii model, the anomalous contribu
The model can be initialized in the configuration as

```julia
conductivity_model = Mitchner(anom_scale = 1.0)
```

In the user-defined configuration


## Custom thermal conductivity models

The procedure for adding additional thermal conductivity models generally follows that of adding anomalous transport models, as described on the [Anomalous transport](@ref) page.
The two main differences is that the model type is of 'ThermalConductivityModel' rather than 'AnomalousTransportModel' and that the function needs to return ``\kappa`` rather than ``\nu an``.
For demonstration pursposes however, let's say we wanted to implement the Braginskii model but without the anomalous collision frequency and a coefficient of 4.66 (singly charged ions only).
We first define the type:

```julia
using HallThruster

struct Braginskii_Classical <: ThermalConductivityModel end
```

With the type defined, we then need to define the function that describes the model. This function should take two arguements, the vector of thermal conductivity values and the solver params.

```julia
function (model::Braginskii_Classical)(κ, params)
    (;B, ne, Tev, νc) = params.cache
    me = HallThruster.me
    e = HallThruster.e

    for i in eachindex(κ)
        ωce = e * B[i] / me
        κ[i] = 4.66 * ne * e * Te * νc[i] / (me * ωce^2)
    end

    return κ
end
```

We could now employ this model by setting `conductivity_model = Braginskii_Classical()' in the user-defined configuration. 
