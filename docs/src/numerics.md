# Numerics

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
See [Fluxes](@ref) for the implemented fluxes, and possible limiters to be used in reconstruction to ensure a total variation diminishing scheme (TVD).

## Time discretization of heavy species

Time integration is handled by [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). Strong stability preserving Runge Kutta schemes are used by default (specifically `SSPRK22`), but the user is free use other schemes provided by DifferentialEquations.jl.

Keep in mind that following von-Neumann stability analysis for explicit schemes the CFL number has to be lower than 1. The CFL number is defined as:

```math
    \sigma = \frac{u_i \Delta t}{\Delta x}
```

Using a maximum ion velocity of 22000 m/s, a domain length of 0.05m and 200 cells, results in ``\Delta t \leq 1.2 \times 10^{-8}`` s. In practice, it needs to be a bit lower in order to handle transients as the solution oscillates. This restriction is valid for the continuity and isothermal euler equations. Information on setting `dt` and selecting the integration scheme can be found in the [Tutorial](@ref).


## Electron energy equation discretization

While the ions can be explicitly solved with ``$\Delta t \sim 10^{-8}``s, the heat condution term in the electron energy equation adds additional constraints which would lower the timestep by about a factor of 10. In order to not further increase the timestepping restrictions and increase computation time, the electron energy equation is solved semi-implicitly in time using a backward Euler or Crank-Nicholson scheme. See [Configuration](@ref) for information on how to select which scheme is used.

The spatial discretization of the electron energy equation uses central finite differences in a manner similar to the potential solver (see below). This, combined with the semi-implicit timestepping, creates a tridiagonal linear system which can be efficiently solved using the Thomas algorithm.


## Potential solver

The second order elliptic differential equation for the potential is discretized using a second order centered difference scheme, with all derivatives appearing approximated in a similar fashion. The values are solved for on the edges of the fluid cells, leading to a staggered grid and avoiding interpolation for the most part. The tridiagonal system ``\underline{\underline{A}} \underline{x} = \underline{b}`` is then solved. The left hand side of the potential equation, see [Physics model](@ref) is discretized as follows, the indexing refers to the fluid discretization: 

```math
    \frac{\partial}{\partial z}\left(\mu_{\perp} n_e \frac{\partial\phi}{\partial z}\right)\bigg\vert^\delta_{_{i + \frac{1}{2}}} \approx \frac{1}{h} \left(\left(\mu_{\perp}\vert^\delta_{i + 1} n_e\vert^\delta_{i + 1} \frac{\partial\phi}{\partial z}\bigg\vert^\delta_{_{i+1}}\right) - \left(\mu_{\perp}\vert^\delta_{i} n_e\vert^\delta_{i} \frac{\partial\phi}{\partial z}\bigg\vert^\delta_{_{i}}\right) \right) + O(h^2)
```

where

```math
    \frac{\partial\phi}{\partial z}\bigg\vert^\delta_{_{i+1}} \approx \frac{1}{h}\left(\phi_{i+\frac{3}{2}} - \phi_{i + \frac{1}{2}}\right) + O(h^2)
```

similarly

```math
    \frac{\partial\phi}{\partial z}\bigg\vert^\delta_{_{i}} \approx \frac{1}{h}\left(\phi_{i + \frac{1}{2}} - \phi_{i-\frac{1}{2}}\right) + O(h^2)
```

results in

```math
    \frac{\partial}{\partial z}\left(\mu_{\perp} n_e \frac{\partial\phi}{\partial z}\right)\bigg\vert^\delta_{_{i}} \approx \frac{1}{h^2} \left(\mu_{\perp}\vert^\delta_{i + 1} n_e\vert^\delta_{i + 1} \phi_{i + \frac{3}{2}} - (\mu_{\perp}\vert^\delta_{i + 1} n_e\vert^\delta_{i + 1} + \mu_{\perp}\vert^\delta_{i} n_e\vert^\delta_{i}) \phi_{i + \frac{1}{2}} + \mu_{\perp}\vert^\delta_{i} n_e\vert^\delta_{i} \phi_{i - \frac{1}{2}}\right)  + O(h^2)
```

The left hand side is incorporated into a ``\mathrm{NxN}`` matrix A, while the RHS is added to the vector b and ``N = n_{cells} + 1``. This results in a tridiagonal matrix that is diagonally dominant, and can therefore be solved using the Thomas algorithm at a computational expense of O(``N``), rather than O(``N^3``) for standard Gaussian elimination. The implementation can be found [here](https://um-pepl.github.io/HallThruster.jl/dev/internals/#HallThruster.solve_potential_edge!-Tuple{Any,%20Any}). 

## Evaluation of derivatives

Some computations require the numerical approximation of derivatives, for example the evaluation of the electron velocity from the equation for electron current using the generalized Ohm's law, see [Physics model](@ref). The derivatives are evaluated to second order using [forward difference](https://um-pepl.github.io/HallThruster.jl/dev/internals/#HallThruster.forward_difference-NTuple{6,%20Any}), [central difference](https://um-pepl.github.io/HallThruster.jl/dev/internals/#HallThruster.central_difference-NTuple{6,%20Any}) or [backward difference](https://um-pepl.github.io/HallThruster.jl/dev/internals/#HallThruster.backward_difference-NTuple{6,%20Any}) depending on the location in the domain. 

## Electron pressure coupled method

The generalized Ohm's law described in [Physics Model](@ref) can be rewritten in the following assuming only one charge state for simplicity.

```math
E = \frac{-u_{e}}{\mu} - \frac{1}{en_i}\nabla (n_i k_B T_e)
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
