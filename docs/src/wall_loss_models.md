# Wall Loss Models

HallThruster.jl allows you to choose from three different wall loss models. They approximate the electron energy lost to the thruster walls in radial direction. As the computational axis of the 1D code is axially in the thruster, the wall loss is not directly resolved by the fluid and applied in each cell as an electron energy loss term. 

## Background

The core of the wall loss models in HallThruster.jl is the abstract type `WallLossModel`. It has three struct children: `NoWallLosses`, `ConstantSheathPotential`, and `WallSheath`. `ConstantSheathPotential` has three fields. A wall `sheath_potential` to be set by the user, and an `inner_loss_coeff` and `outer_loss_coeff` which allow to scale the energy loss inside vs. outside the thruster channel. `WallSheath` has one field `material`, which is of struct `WallMaterial`.

## Provided wall loss models

HallThruster.jl provides three models out of the box. These are

| Model                   | Supported species                                            | Description                                                  |
| ----------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| `NoWallLosses` | `Any`                                                      | Ignores electron energy losses to the walls. May cause numerical issues.  |
| `ConstantSheathPotential`         | `Any` | Employs a simple sheath energy loss model with constant sheath potential, based on the electron Boltzmann equation for electron density in the sheath as a function of electron temperature. Uses constants to scale losses inside and outside the thruster. See also [JP Boeuf, *Low frequency oscillations in a stationary plasma thruster*, Journal of Applied Physics 84, 3541, 1998](https://aip.scitation.org/doi/10.1063/1.368529) and [Landmark study](https://www.landmark-plasma.com/test-case-3)|
| `WallSheath`         | `Xenon`                                                      | Conceputally similar loss model as `ConstantSheathPotential`, but evaluates constants and sheath potential given in the previously mentioned using approximations. The [effective collision frequency](https://um-pepl.github.io/HallThruster.jl/dev/internals/#HallThruster.effective_loss_frequency-Tuple{Any}) follows from a half-maxwellian approximation.  The [sheath potential](https://um-pepl.github.io/HallThruster.jl/dev/internals/#HallThruster.compute_wall_sheath_potential-Tuple{Any,%20Any,%20Any})  is evaluated as a function of [secondary electron emission](https://um-pepl.github.io/HallThruster.jl/dev/internals/#HallThruster.SEE_yield-Tuple{HallThruster.WallMaterial,%20Any}), which is a function of wall material, and electron temperature. This model is based on [Hobbs and Wesson, *Heat flow through a Langmuir sheath in the presence of electron emission*, Plasma Physics 9 85, 1967](https://iopscience.iop.org/article/10.1088/0032-1028/9/1/410) and described in [Goebel and Katz, *Fundamentals of Electric Propulsion*, JPL Science and Technology series, 2008](https://descanso.jpl.nasa.gov/SciTechBook/series1/Goebel__cmprsd_opt.pdf)

The density relation in the 'WallSheath' model is based upon the electron [`Boltzmann relation`](@ref). 
