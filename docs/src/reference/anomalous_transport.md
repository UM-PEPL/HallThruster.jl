# Anomalous transport

HallThruster has a few anomalous transport models built in and allows users to define their own. This page describes these models and the process by which algebraic and multi-equation transport
models can be added by the user.

## Built-in Models

```@meta
CurrentModule = HallThruster
```

```@docs
NoAnom
Bohm
TwoZoneBohm
GaussianBohm
MultiLogBohm
LogisticPressureShift
```

## The `AnomalousTransportModel` interface

Users defining their own transport model will need first define a model that subtypes `AnomalousTransportModel`.
They will then need to provide definitions for the following methods

```@docs
num_anom_variables
```
Additionally, they will need to define a function that allows the model to be called with the following arguments.

```julia
# replace <:AnomalousTransportModel with ::MyModel, where MyModel is the name of your model.
(model<:AnomalousTransportModel)(nu_an, params, config)
```

This function should operate in-place on nu_an.
See [Adding an anomalous transport model](@ref) for a guide on implementing new models.
