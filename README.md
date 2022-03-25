[![HallThruster.jl](https://raw.githubusercontent.com/archermarx/HallThruster.jl/main/docs/src/assets/banner.svg)](https://archermarx.github.io/HallThruster.jl/dev)


| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://archermarx.github.io/HallThruster.jl/dev) | [![CI](https://github.com/archermarx/HallThruster.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/archermarx/HallThruster.jl/actions/workflows/ci.yml) [![codecov](https://codecov.io/gh/archermarx/HallThruster.jl/branch/main/graph/badge.svg?token=cEoGN49eZp)](https://codecov.io/gh/archermarx/HallThruster.jl)|


HallThruster.jl is an open-source, 1D fluid Hall thruster code written in Julia. It was initially developed by the University of Michigan's [Plasmadynamics and Electric Propulsion Laboratory](https://pepl.engin.umich.edu) and licensed under the MIT license.

## Installation

To install HallThruster.jl, you must first install Julia 1.7 or above from the [official Julia site](https://julialang.org/downloads/), or by using [juliaup](https://github.com/JuliaLang/juliaup). We recommend using the latest Julia release when possible. Once installed, launch Julia and type

```julia
julia> ]add https://github.com/archermarx/HallThruster.jl
```

This will install HallThruster.jl using Julia's package manager. For details on setting up and running Hall thruster simulations, see [the official documentation](https://archermarx.github.io/HallThruster.jl/dev).

## Physics model

HallThruster.jl solves the quasineutral plasma equations of motion for a Hall Thruster along the thruster's channel centerline (the z-axis). We solve seperate models for neutral particles, ions, and electrons. Neutrals are assumed to have (user-configurable) constant velocity and temperature and are tracked by a single continuity equation. Ions are assumed isothermal and unmagnetized. Multiple ion species with different charge states are supported, and each is tracked by a continuity equation and a momentum equation. We employ the drift-diffusion approximation for electrons, which reduces the electron momentum equation to a generalized Ohm's law. Charge conservation is then used to solve for the electrostatic potential. The electron temperature is determined by solving an equation for the conservation of electron internal energy. 

For more details as well as several references, see the documentation.
