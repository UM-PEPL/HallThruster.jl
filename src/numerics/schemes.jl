"""
    $(TYPEDEF)
Define a hyperbolic scheme for the heavy species solve.
Consists of a `FluxFunction`, a `SlopeLimiter`, and a `Bool`,
the latter indicating whether gradient-reconstruction should be employed to obtain
second-order accuracy
"""
@kwdef struct HyperbolicScheme{F, L}
    flux_function::FluxFunction{F} = rusanov
    limiter::SlopeLimiter{L} = van_leer
    reconstruct::Bool = true
end
