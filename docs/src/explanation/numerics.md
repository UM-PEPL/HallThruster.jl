# Numerics

HallThruster.jl uses the finite volume method and the Rusanov flux for calculating ion and neutral fluid fluxes.
Second-order gradient reconstruction is available using the van Leer limiter.
Time-marching for the heavy species is handled using a second order strong-stability preserving Runge-Kutta scheme (SSPRK22).
The left hand side of the electron energy equation is integrated implicitly using a Crank Nicolson Adams Bashforth (CNAB) scheme.
This enables larger timesteps due to the severe restrictions due to the electron heat flux.

## Spatial discretization for heavy species

Neutrals and ions are considered heavy species (compared to the much lighther electrons). HallThruster.jl uses the finite volume method (FVM). FVM has the advantage that it is by definition conservative, which is a useful property when solving hyperbolic conservation laws such as the Euler equations.
Currently, only the continuity equation is solved for the neutrals and the isothermal Euler equations for the ion species.

```math
\int_{i-\frac{1}{2}}^{i+\frac{1}{2}} \frac{\partial n_n}{\partial t} \,dz + \int_{i-\frac{1}{2}}^{i+\frac{1}{2}} \frac{\partial n_n u_n}{\partial z} \,dz = \int_{i-\frac{1}{2}}^{i+\frac{1}{2}} \dot{n}_n \, dz
```

The ``n_n u_n`` can be replaced by a generic flux term ``F(z)`` and generalized to any advection like equation. Treatment of the source term is described in [Collisions and Reactions](@ref). Integration results in

```math
h\frac{\partial n_n}{\partial t} + \left(F_{_{i+\frac{1}{2}}} - F_{_{i-\frac{1}{2}}}\right) = h \dot{n}_n
```

## Time discretization of heavy species

We employ a strong stability preserving second-order Runge Kutta scheme (`SSPRK22`) for timestepping

The user has the option to supply a fixed timestep, or a CFL number. In the former case, the user will need to select a timestep that obeys the CFL condition, defined as

```math
    \sigma = \frac{u_i \Delta t}{\Delta x} < 1
```

Using a maximum ion velocity of 22000 m/s, a domain length of 0.05m and 200 cells, results in ``\Delta t \leq 1.2 \times 10^{-8}`` s.
In practice, it needs to be a bit lower in order to handle transients as the solution oscillates. This restriction is valid for the continuity and isothermal euler equations.

In most cases, it is better to let HallThruster.jl handle timestepping automatically using its adaptive timestepping option. If adaptive timestepping is enabled, the user-defined timestep is ignored in favor of a timestep based on the minimum of three conditions and a user-supplied CFL number. Mathematically the timstep is choosen as:

```math
    \Delta t = min(\sigma \frac{\Delta x}{max(u_i + a_i, u_i - a_i)}, \sigma \frac{\dot{n}_i}{n_i}, \sqrt{\frac{\sigma m_i \Delta x}{q_i E}})
```

Where ``a_i`` is the ion sound speed. Physically, these three conditions represent timestep limits imposed by the flux, ionization, and electrostatic acceleration. Keep in mind that due to stability limits imposed by the ionization condition, the CFL number cannot be higher than 0.799 to remain stable. This limit will be imposed by HallThruster.jl if the user-defined value is too high.

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
