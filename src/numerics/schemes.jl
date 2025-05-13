"""
    $(TYPEDEF)
Define a hyperbolic scheme for the heavy species solve.
Consists of a `FluxFunction`, a `SlopeLimiter`, and a `Bool`,
the latter indicating whether gradient-reconstruction should be employed to obtain
second-order accuracy
"""
@kwdef struct HyperbolicScheme{L}
    limiter::SlopeLimiter{L} = van_leer
    flux_function::Symbol = :rusanov
    reconstruct::Bool = true
end
