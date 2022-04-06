# Fluxes and Numerics

As described in [Configuration](@ref) and [Initialization](@ref) different flux options are available in `HyperbolicScheme`. Timemarching for the heavy species is handled by [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). The left hand side of the electron energy equation is integrated implicitly using a Crank Nicolson Adams Bashforth (CNAB) scheme. This enables larger timessteps due to the severe restrictions due to the electron heat flux. 

## Spatial discretization for heavy species

Neutrals and ions are considered heavy species (compared to electrons). HallThruster.jl uses the finite volume method (FVM). FVM has the advantage that it is by definition conservative, which is a useful property when solving hyperbolic conservation laws such as the Euler equations. Currently, only the continuity equation is solved for the neutrals and the isothermal Euler equations for the ion species. Possibly, the full Euler equations will be added in the future, its implementation has been verified using the Sod Shock tube. The following provides and example of the control volume approach applied to the continuity equation. 

```math
\int_{i-\frac{1}{2}}^{i+\frac{1}{2}} \frac{\partial n_n}{\partial t} \,dz + \int_{i-\frac{1}{2}}^{i+\frac{1}{2}} \frac{\partial n_n u_n}{\partial z} \,dz = \int_{i-\frac{1}{2}}^{i+\frac{1}{2}} \dot{n}_n \, dz
```

The ``n_n u_n`` can be replaced by a generic flux term ``F(z)`` and generalized to any advection like equation. Treatment of the source term is described in [Collisions and Reactions](@ref). Integration results in

```math
h\frac{\partial n_n}{\partial t} + \left(F_{_{i+\frac{1}{2}}} - F_{_{i-\frac{1}{2}}}\right) = h \dot{n}_n
```


The fluxes ``F_{_{i+\frac{1}{2}}}`` and ``F_{_{i-\frac{1}{2}}}`` are reconstructed at the cell interfaces, and for this flux reconstruction multiple options are available. These are set using the object `HyperbolicScheme` consisting of fields `flux`, `limiter`, `reconstruct` and `WENO`.

Three different flux approximations are available. 

| Flux                   |  Description                                                  |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| `upwind`                                                    | Simple first order accurate flux approximation, that as a results does not distinguish between cell centered and cell average values and adapts reconstruction according to sign of advection velocity. Very diffusive. No Riemann solver or approximation.   |
| `HLLE`       | Approximate Riemann solver. The Harten-Lax-van Leer-Einfeldt scheme approximates a Riemann problem with three constant states. see reference. The scheme is positively-conservative if stability bounds for maximum and minimum wavespeeds are met, which makes it useful in its application with HallThruster.jl. First order accurate in space. [B. Einfeldt. *On godunov-type methods for gas dynamics.* Journal of Computational Physics, 25:294-318, 1988.](https://epubs.siam.org/doi/10.1137/0725021) |
| `rusanov`                                                     | Approximate Riemann solver. Also known as the local Lax-Friedrich flux. Has slighlty modified choice of wave speeds. Adds viscosity to a centered flux. More diffusive than HLLE. [Chi-Wang Shu, *Lecture Notes: Numerical Methods for Hyperbolic Conservation Laws (AM257)*](https://mathema.tician.de/dl/academic/notes/257/257.pdf)


These flux approximations are all first order accurate in space (piecewise constant recontruction), but can be extended to piecewise linear reconstruction within a cell. To satisfy stability bounds and keep the scheme total variation diminishing (TVD), it has to be coupled with a limiter. Many limiters have been proposed, the ones implemented in HallThruster.jl are the following: `koren`, `minmod`, `osher`, `superbee`, `van_albada`, `van_leer`. If the field `reconstruction` is set to `true`, the selected limiter will be used. 


`WENO` refers to a 5th order weighted essentially non-oscillatory that will be added but not yet functional and should be set to `false`. Depending on the smoothness of the solution, it is capable of 5th order spacial accuracy HLLE flux as a building block. [Chi-Wang Shu. *High-order finite difference and finite volume weno schemes and discontinuous galerkin
methods for cfd.* International Journal of Computational Fluid Dynamics, 17(2):107–118, 2003](https://www.tandfonline.com/doi/abs/10.1080/1061856031000104851) [Dinshaw S. Balsara and Chi-Wang Shu. Monotonicity preserving weighted essentially non-oscillatory schemes with increasingly high order of accuracy. Journal of Computational Physics, 160(2):405–452, 2000.](https://www.sciencedirect.com/science/article/pii/S002199910096443X)


## Time discretization heavy species

Time integration is handled by [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). Strong stability preserving Runge Kutta schemes are used by default (specifically `SSPRK22`), but the user is free use other schemes provided by DifferentialEquations.jl.

Keep in mind that following von-Neumann stability analysis for explicit schemes the CFL number has to be lower than 1. The CFL number is defined as:

```math
    \sigma = \frac{u_i \Delta t}{\Delta x}
```

Using a maximum ion velocity of 22000 m/s, a domain length of 0.05m and 200 cells, results in ``\Delta t \leq 1.2e-8`` s. This is valid for the continuity and isothermal euler equations, due to the heat flux term in the electron energy equation an additional constraint is added. In order to not further increase the timestepping restrictions and increase computational complexity, the electron energy equation is solved semi-implicitly. Setting `dt` and selecting the integration scheme is shown in [Initialization](@ref). 


## Electron energy equation discretization


## Potential solver

The second order elliptic differential equation for the potential is solved using a finite difference approach. The potential is solved for on the edges of the fluid cells, which avoids the necessity of interpolation. A first order central difference scheme is used for the first derivative, and a second order central difference scheme for the second derivative. The discretization can be formulated as a tridiagonal system ``\underline{\underline{A}} \underline{x} = \underline{b}``. This is solved using the Thomas' algorithm, which reduces computational complexity to the order O(N). More information can be found in the function [definition](https://um-pepl.github.io/HallThruster.jl/dev/internals/#HallThruster.solve_potential_edge!-Tuple{Any,%20Any}). 

## Evaluation of derivatives

Some computations require the numerical approximation of derivatives, for example the evaluation of the electron velocity from the equation for electron current using the generalized Ohm's law, see [Physics model](@ref). The derivatives are evaluated to second order using [forward difference](https://um-pepl.github.io/HallThruster.jl/dev/internals/#HallThruster.forward_difference-NTuple{6,%20Any}), [central difference](https://um-pepl.github.io/HallThruster.jl/dev/internals/#HallThruster.central_difference-NTuple{6,%20Any}) or [backward difference](https://um-pepl.github.io/HallThruster.jl/dev/internals/#HallThruster.backward_difference-NTuple{6,%20Any}) depending on the location in the domain. 

## Electron pressure coupled method

The generalized Ohm's law described in [Physics Model](@ref) can be rewritten in the following assuming only one charge state for simplicity.

```math
E &= \frac{-u_{e}}{\mu} - \frac{1}{en_i}\nabla (n_i k_B T_e)
```

Substituting above expression into the ion momentum equation results in 

```math
    \frac{\partial}{\partial t} (n_{i} u_{i}) + \frac{\partial}{\partial z} (n_{i} u_{i}^2 + \frac{p_{i}}{m_i}) = \frac{e}{m_i}\left(n_{i} \frac{-u_{e}}{\mu} - \frac{1}{en_i}\nabla (n_i k_B T_e)\right)
```

and can be rearranged into 

```math
    \frac{\partial}{\partial t} (n_{i} u_{i}) + \frac{\partial}{\partial z} (n_{i} u_{i}^2 + \frac{p_{i} + p_{e}}{m_i}) = \frac{-e n_{i} u_{e}}{m_i \mu}
```

This effectively increases the speed of sound in the ion momentum equation from the ion thermal speed to the ion acoustic speed, which decreases the Mach range covered in the domain and provides a direct feedback from the electron pressure to the ion momentum, thereby suppressing numerical oscillations. It is recommended to set `electron_pressure_coupled` to `1.0` in [Configuration](@ref). See [K. Hara, *Non-oscillatory quasineutral fluid model of cross-field discharge plasmas*, Physics of Plasmas 25, 123508, 2018](https://aip.scitation.org/doi/pdf/10.1063/1.5055750). 
