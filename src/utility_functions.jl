"""
    smooth_max(x, y, k=10)
Computes a smooth approximation to max(x, y)
"""
function smooth_max(x,y,k = 10)
    (x*exp(k*x) + y*exp(k*y)) / (exp(k*x) + exp(k*y))
end

"""
    smooth_min(x, y, k=10)
Compute a smooth approximation to min(x, y)
"""
smooth_min(x,y,k=10) = smooth_max(x, y, -k)

"""
    smooth_if_gt(x, cutoff, v1, v2, k=10)
Computes an analytic approximation to x < cutoff ? v1 : v2
"""
smooth_if(x, cutoff, v1, v2, k = 10) = 0.5*((v2-v1)*tanh(k*(x-cutoff)) + v1+v2)

"""
    mutable struct EnergyOVS
Enables setting mu, ue, Tev and ne to certain values to very electron energy equation
"""
mutable struct EnergyOVS{F1, F2}
    active ::Int64
    μ::Union{Float64, Nothing}
    ue::Union{Float64, Nothing}
    Tev::F1
    ne::F2
end

"""
    mutable struct Verification
is passed to params to identify if OVS is active.
"""

mutable struct Verification{F1, F2}
    potential ::Int64
    fluid ::Int64
    energy ::EnergyOVS{F1, F2}
end

@inline function uneven_forward_coeffs(z0, z1, z2)
    h1 = z1 - z0
    h2 = z2 - z1
    return (
        -(2h1 + h2)/h1/(h1 + h2),
        (h1 + h2)/(h1*h2),
        - h1/h2/(h1 + h2)
    )
end

@inline function uneven_central_coeffs(z0, z1, z2)
    h1 = z1 - z0
    h2 = z2 - z1
    return (
        -h2/h1/(h1+h2),
        -(h1-h2)/(h1*h2),
        h1/h2/(h1+h2)
    )
end

@inline function uneven_backward_coeffs(z0, z1, z2)
    h1 = z1 - z0
    h2 = z2 - z1
    return (
        h2/h1/(h1 + h2),
        -(h1 + h2)/(h1*h2),
        (h1+2h2)/h2/(h1+h2),
    )
end

@inline function uneven_second_deriv_coeffs(z0, z1, z2)
    h1 = z1 - z0
    h2 = z2 - z1
    common = 2 / (h1 * h2 * (h1 + h2))
    return (
        h2 * common,
        -(h1 + h2) * common,
        h1 * common,
    )
end

@inline function uneven_forward_diff(f0, f1, f2, z0, z1, z2)
    c0, c1, c2 = uneven_forward_coeffs(z0, z1, z2)
    return c0 * f0 + c1 * f1 + c2 * f2
end

@inline function uneven_central_diff(f0, f1, f2, z0, z1, z2)
    c0, c1, c2 = uneven_central_coeffs(z0, z1, z2)
    return c0 * f0 + c1 * f1 + c2 * f2
end

@inline function uneven_backward_diff(f0, f1, f2, z0, z1, z2)
    c0, c1, c2 = uneven_backward_coeffs(z0, z1, z2)
    return c0 * f0 + c1 * f1 + c2 * f2
end

@inline function uneven_second_deriv(f0, f1, f2, z0, z1, z2)
    c0, c1, c2 = uneven_second_deriv_coeffs(z0, z1, z2)
    return c0 * f0 + c1 * f1 + c2 * f2
end

@inline function diff(f0, f1, z0, z1)
    return (f1 - f0) / (z1 - z0)
end
