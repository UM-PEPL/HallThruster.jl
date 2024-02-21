# Electron Thermal Conductivity

HallThruster has a few electron thermal conductivity models built in and allows users to define their own. This page describes these models and the process by which algebraic transport
models can be added by the user.

## Built-in Models

### `Mitchner()`

HallThruster.jl's default electron thermal conductivity option. This model employs the form described in Mitchner & Kruger "Partially Ionized Gases" (1973). The thermal conductivity takes the form:

```math
\begin{aligned}
    \kappa_{e\perp} & \approx \frac{1}{1+\Omega_H^2} \frac{2.4}{1+ \frac{\nu_{ie}}{\sqrt{2}\nu}} \frac{neTe (eV)}{\nu m_e} \\
    \nu = \nu_{ei} + \nu_{en} + \nu_{AN} \\
\end{aligned}
```

The model is initialized by default but can be explicitly enabled by including:

```julia
conductivity_model = Mitchner()
```

In the user-defined configuration

### `Braginskii()`

This form of the thermal conducitivity follows the result of S. I. Braginskii, inReviews of Plasma Physics, edited byM. A. Leontovich (Consultants Bureau, New York, 1965), Vol. 1,p. 205.:

\begin{aligned}
    \kappa_{e\perp} & \approx C \frac{neTe (eV) \nu}{m_e \omega_{ce}^2} \\
    \nu = \nu_{ei} + \nu_{en} + \nu_{AN} \\
\end{aligned}

Where C is a constant that is based on the value of the effective charge for multiple charge states and Table 1 of the Braginskii reference. A `Braginski` model can be initialized in the user-defined configuration as:

```julia
conductivity_model = Braginskii()
```

## Custom thermal conductivity models

The procedure for adding additional thermal conductivity models generally follows that of adding anomalous transport models, as described on the [Anomalous Transport](@ref) page. The two main differences is that the model type is of 'ThermalConductivityModel' rather than 'AnomalousTransportModel' and that the function needs to return ``\kappa`` rather than ``\nu an``. For demonstration pursposes however, let's say we wanted to implement the Braginskii model but without the anomalous collision frequency and a coefficient of 4.66 (singly charged ions only). We first define the type:

```julia
using HallThruster

struct Braginskii_Classical <: ThermalConductivityModel end
```

With the type defined, we then need to define the function that describes the model. This function should take two arguements, the vector of thermal conductivity values and the solver params.

```julia
function (model::Braginskii_Classical)(κ, params)

    me = HallThruster.me
    e = HallThruster.e
    B = params.cahce.B
    ne = params.cache.ne
    Te = params.cahce.Tev

    for i in eachindex(κ)

        ωce = e * B[i] / me

        κ[i] = 4.66 * ne * e * Te * params.cache.νc[i] / (me * ωce^2)
    end

    return κ
end
```

We could now employ this model by setting `conductivity_model = Braginskii_Classical()' in the user-defined configuration. 
