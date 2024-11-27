![HallThruster.jl Logo](./docs/src/assets/banner_light.svg#gh-dark-mode-only)
![HallThruster.jl Logo](./docs/src/assets/banner.svg#gh-light-mode-only)

| **Documentation** | **Build Status**| **Paper**| **Repository DOI**|
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://UM-PEPL.github.io/HallThruster.jl/dev) | [![CI](https://github.com/UM-PEPL/HallThruster.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/UM-PEPL/HallThruster.jl/actions/workflows/ci.yml) [![codecov](https://codecov.io/gh/UM-PEPL/HallThruster.jl/branch/main/graph/badge.svg?token=cEoGN49eZp)](https://codecov.io/gh/UM-PEPL/HallThruster.jl)| [![status](https://joss.theoj.org/papers/ce9cb7aa54df10d69ed248912e584f53/status.svg)](https://joss.theoj.org/papers/ce9cb7aa54df10d69ed248912e584f53) | [![DOI](https://zenodo.org/badge/394711445.svg)](https://zenodo.org/badge/latestdoi/394711445) |


HallThruster.jl is an open-source, 1D fluid Hall thruster code written in Julia. It is developed by [Thomas Marks](https://thomasmarks.space), [Paul Schedler](https://www.linkedin.com/in/paul-schedler-1b3b6b171/) and Declan Brick at the University of Michigan's [Plasmadynamics and Electric Propulsion Laboratory](https://pepl.engin.umich.edu) and is licensed under the MIT license.

## Installation

To install HallThruster.jl, you must first install Julia 1.7 or above from the [official Julia site](https://julialang.org/downloads/), or by using [juliaup](https://github.com/JuliaLang/juliaup). We recommend using the latest Julia release when possible. Once installed, launch Julia and type `]` to enter the Pkg REPL. To install HallThruster.jl type

```julia
(@v1.10) pkg> add https://github.com/UM-PEPL/HallThruster.jl
```

This will install HallThruster.jl using Julia's package manager. For details on setting up and running Hall thruster simulations, see [the official documentation](https://UM-PEPL.github.io/HallThruster.jl/dev). A Tutorial is available [here](https://nbviewer.org/github/UM-PEPL/HallThruster.jl/blob/main/HallThrusterTutorial.ipynb).

## Contribution

Users are welcome to suggest and implement features for the code, as well as report bugs or numerical issues they encounter. Please feel free to [open an issue on this repository](https://github.com/UM-PEPL/HallThruster.jl/issues/new) describing your desired change/bug-fix. Pull requests are also welcome!

## Citation

If you use this code in your work, please cite our [publication in the Journal of Open Source Software](https://joss.theoj.org/papers/10.21105/joss.04672):

```
@article{Marks2023, 
  doi = {10.21105/joss.04672},
  url = {https://doi.org/10.21105/joss.04672},
  year = {2023},
  publisher = {The Open Journal},
  volume = {8}, number = {86}, pages = {4672},
  author = {Thomas Marks and Paul Schedler and Benjamin Jorns},
  title = {HallThruster.jl: a Julia package for 1D Hall thruster discharge simulation},
  journal = {Journal of Open Source Software} 
} 
```

 This BibTEX entry can also be found in [CITATION.bib](https://github.com/UM-PEPL/HallThruster.jl/blob/main/CITATION.bib).
