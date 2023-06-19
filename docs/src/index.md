```@meta
CurrentModule = HallThruster
DocTestSetup = quote
    using HallThruster
end
```

# HallThruster.jl

HallThruster.jl is an open-source, 1D fluid Hall thruster code written in Julia. It was initially developed by the University of Michigan's [Plasmadynamics and Electric Propulsion Laboratory](https://pepl.engin.umich.edu) and is licensed under the MIT license. 

## Installation

To install HallThruster.jl, you must first install Julia 1.7 or above from the [official Julia site](https://julialang.org/downloads/), or by using [juliaup](https://github.com/JuliaLang/juliaup). We recommend using the latest Julia release when possible. Once installed, launch Julia and type

```julia
julia> ]add https://github.com/UM-PEPL/HallThruster.jl
```

This will install HallThruster.jl using Julia's package manager. For details on setting up and running Hall thruster simulations, see [the official documentation](https://UM-PEPL.github.io/HallThruster.jl/dev). A Tutorial is available [here](https://nbviewer.org/github/UM-PEPL/HallThruster.jl/blob/main/HallThrusterTutorial.ipynb).

## Contribution

Users are welcome to suggest and implement features for the code, as well as report bugs or numerical issues they encounter. Please feel free to [open an issue](https://github.com/UM-PEPL/HallThruster.jl/issues/new) on this repository describing your desired change or bug-fix. Pull requests are also welcome!