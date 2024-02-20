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

HallThruster.jl also has an adaptive timestepping option. If adaptive timestepping is enabled, the user-defined timestep is ignored in favor of a timestep based on the minimum of three conditions and the user-defined CFL number. Mathematically the timstep is choosen as:

```math
    \Delta t = min(\sigma \frac{\Delta x}{max(u_i + a_i, u_i - a_i)}, 0.799*\sigma \frac{\dot{n}_i}{n_i}, \sqrt{\frac{\sigma m_i \Delta x}{q_i E}})
```

Where ``a_i`` is the ion sound speed. Physically, these three conditions represent timestep limits imposed by the flux, ionization, and electrostatic acceleration. The 0.799 factor in front of the ionization term comes from empirical testing, so it may be worth trying a lower CFL if the adaptive timestepping is unstable.  


## Electron energy equation discretization

While the ions can be explicitly solved with ``\Delta t \sim 10^{-8}``s, the heat condution term in the electron energy equation adds additional constraints which would lower the timestep by about a factor of 10. In order to not further increase the timestepping restrictions and increase computation time, the electron energy equation is solved semi-implicitly in time using a backward Euler or Crank-Nicholson scheme. See [Configuration](@ref) for information on how to select which scheme is used.

The spatial discretization of the electron energy equation uses central finite differences in a manner similar to the potential solver (see below). This, combined with the semi-implicit timestepping, creates a tridiagonal linear system which can be efficiently solved using the Thomas algorithm.

If adaptive timestepping is enabled, the timestep used for the explicit terms is limited by: 

```math
    \Delta t = \abs{\frac{3\sigma n_eT_e}{W_{loss} + S_{coll}}}
```

As this timestep may be significantly smaller than the timestep used for the heavy species, the electron energy equation is updated using a series of sub-steps with the ``\Delta t `` enforced by the equation above until a total ``\Delta t `` equal to that set by the heavy species is reached. 

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
