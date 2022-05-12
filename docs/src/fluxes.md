# Fluxes

HallThruster.jl uses the Finite Volume method, and as such the face values of the fluxes need to be reconstructed. See [Numerics](@ref) for more information. 


The fluxes ``F_{_{i+\frac{1}{2}}}`` and ``F_{_{i-\frac{1}{2}}}`` are reconstructed at the cell interfaces, and for this flux reconstruction multiple options are available. These are set using the object `HyperbolicScheme` consisting of fields `flux`, `limiter`, `reconstruct` and `WENO`.

Three different flux approximations are available. 

| Flux                   |  Description                                                  |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| `upwind`                                                    | Simple first order accurate flux approximation, that as a results does not distinguish between cell centered and cell average values and adapts reconstruction according to sign of advection velocity. Very diffusive. No Riemann solver or approximation.   |
| `HLLE`       | Approximate Riemann solver. The Harten-Lax-van Leer-Einfeldt scheme approximates a Riemann problem with three constant states. see reference. The scheme is positively-conservative if stability bounds for maximum and minimum wavespeeds are met, which makes it useful in its application with HallThruster.jl. First order accurate in space. [B. Einfeldt. *On godunov-type methods for gas dynamics.* Journal of Computational Physics, 25:294-318, 1988.](https://epubs.siam.org/doi/10.1137/0725021) |
| `rusanov`                                                     | Approximate Riemann solver. Also known as the local Lax-Friedrich flux. Has slighlty modified choice of wave speeds. Adds viscosity to a centered flux. More diffusive than HLLE. [Chi-Wang Shu, *Lecture Notes: Numerical Methods for Hyperbolic Conservation Laws (AM257)*](https://mathema.tician.de/dl/academic/notes/257/257.pdf)


These flux approximations are all first order accurate in space (piecewise constant recontruction), but can be extended to piecewise linear reconstruction within a cell. To satisfy stability bounds and keep the scheme total variation diminishing (TVD), it has to be coupled with a limiter. Many limiters have been proposed, the ones implemented in HallThruster.jl are the following: `koren`, `minmod`, `osher`, `van_albada`, `van_leer`. If the field `reconstruction` is set to `true`, the selected limiter will be used.

`WENO` refers to a 5th order weighted essentially non-oscillatory that will be added but not yet functional and should be set to `false`. Depending on the smoothness of the solution, it is capable of 5th order spacial accuracy HLLE flux as a building block. [Chi-Wang Shu. *High-order finite difference and finite volume weno schemes and discontinuous galerkin
methods for cfd.* International Journal of Computational Fluid Dynamics, 17(2):107–118, 2003](https://www.tandfonline.com/doi/abs/10.1080/1061856031000104851) [Dinshaw S. Balsara and Chi-Wang Shu. Monotonicity preserving weighted essentially non-oscillatory schemes with increasingly high order of accuracy. Journal of Computational Physics, 160(2):405–452, 2000.](https://www.sciencedirect.com/science/article/pii/S002199910096443X)
