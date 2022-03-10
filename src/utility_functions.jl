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


"""
    forward_diff_coeffs(x0, x1, x2)
Generate finite difference coefficients for a forward first derivative approximation  at the point x0
on a three-point stencil at points x0, x1, and x2

```jldoctest;setup = :(using HallThruster: forward_diff_coeffs)
julia> forward_diff_coeffs(1.0, 2.0, 3.0)
(-1.5, 2.0, -0.5)
julia> forward_diff_coeffs(0//1, 1//2, 3//2)
(-8//3, 3//1, -1//3)
```
"""
@inline function forward_diff_coeffs(x0, x1, x2)
    h1 = x1 - x0
    h2 = x2 - x1
    return (
        -(2h1 + h2)/h1/(h1 + h2),
        (h1 + h2)/(h1*h2),
        - h1/h2/(h1 + h2)
    )
end

"""
    central_diff_coeffs(x0, x1, x2)
Generate finite difference coefficients for a central first derivative approximation at the point x1
on a three-point stencil at points x0, x1, and x2

```jldoctest;setup = :(using HallThruster: central_diff_coeffs)
julia> central_diff_coeffs(-1//1, 0//1, 1//1)
(-1//2, 0//1, 1//2)
julia> central_diff_coeffs(-1//2, 0//1, 1//1)
(-4//3, 1//1, 1//3)
```
"""
@inline function central_diff_coeffs(x0, x1, x2)
    h1 = x1 - x0
    h2 = x2 - x1
    return (
        -h2/h1/(h1+h2),
        -(h1-h2)/(h1*h2),
        h1/h2/(h1+h2)
    )
end

"""
    backward_diff_coeffs(x0, x1, x2)
Generate finite difference coefficients for a backward first derivative approximation at the point x2
on a three-point stencil at points x0, x1, and x2

```jldoctest;setup = :(using HallThruster: backward_diff_coeffs)
julia> backward_diff_coeffs(-2//1, -1//1, 0//1)
(1//2, -2//1, 3//2)
julia> backward_diff_coeffs(-3//2, -1//1, 0//1)
(4//3, -3//1, 5//3)
```
"""
@inline function backward_diff_coeffs(x0, x1, x2)
    h1 = x1 - x0
    h2 = x2 - x1
    return (
        h2/h1/(h1 + h2),
        -(h1 + h2)/(h1*h2),
        (h1+2h2)/h2/(h1+h2),
    )
end

"""
    second_deriv_coeffs(x0, x1, x2)
Generate finite difference coefficients for a central second derivative approximation at the point x1
on a three-point stencil at points x0, x1, and x2

```jldoctest;setup = :(using HallThruster: second_deriv_coeffs)
julia> second_deriv_coeffs(-2//1, 0//1, 2//1)
(1//4, -1//2, 1//4)
julia> second_deriv_coeffs(-1//2, 0//1, 1//1)
(8//3, -12//3, 4//3)
```
"""
@inline function second_deriv_coeffs(x0, x1, x2)
    h1 = x1 - x0
    h2 = x2 - x1
    common = 2 / (h1 * h2 * (h1 + h2))
    return (
        h2 * common,
        -(h1 + h2) * common,
        h1 * common,
    )
end

"""
    forward_difference(f0, f1, f2, x0, x1, x2)
Given three points x0, x1, and x2, and the function values at those points, f0, f1, f2,
compute the second-order approximation of the derivative at x0

```jldoctest;setup = :(using HallThruster: forward_difference)
julia> f(x) = x^4; forward_difference(f(2), f(2 + 2*eps(Float64)), f(2 + 2*eps(Float64));
```
"""
@inline function forward_difference(f0, f1, f2, x0, x1, x2)
    c0, c1, c2 = forward_diff_coeffs(x0, x1, x2)
    return c0 * f0 + c1 * f1 + c2 * f2
end

@inline function central_difference(f0, f1, f2, x0, x1, x2)
    c0, c1, c2 = central_diff_coeffs(x0, x1, x2)
    return c0 * f0 + c1 * f1 + c2 * f2
end

@inline function backward_difference(f0, f1, f2, x0, x1, x2)
    c0, c1, c2 = backward_diff_coeffs(x0, x1, x2)
    return c0 * f0 + c1 * f1 + c2 * f2
end

@inline function second_deriv_central_diff(f0, f1, f2, x0, x1, x2)
    c0, c1, c2 = second_deriv_coeffs(x0, x1, x2)
    return c0 * f0 + c1 * f1 + c2 * f2
end

@inline function diff(f0, f1, x0, x1)
    return (f1 - f0) / (x1 - x0)
end
