"""
    smooth_if(x, cutoff, v1, v2, k=10)
Computes an analytic approximation to x < cutoff ? v1 : v2

```jldoctest;setup = :(using HallThruster: smooth_if)
v1 = 1
v2 = 2
cutoff = 0
smooth_if(10.0, cutoff, v1, v2) ≈ v2

# output

true

````
"""
smooth_if(x, cutoff, v1, v2, k = 10) = 0.5*((v2-v1)*tanh(k*(x-cutoff)) + v1+v2)

function smooth_transition(x, x_exh, L_trans, a_in, a_out)
    return a_in + (a_out - a_in) * (x - x_exh)^2 / ((x - x_exh)^2 + L_trans^2)
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
(8//3, -4//1, 4//3)
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
f(x) = x^4
x0, x1, x2 = 2.0, 2.000001, 2.000002
forward_difference(f(x0), f(x1), f(x2), x0, x1, x2)

# output

31.999999998137355
```
"""
@inline function forward_difference(f0, f1, f2, x0, x1, x2)
    c0, c1, c2 = forward_diff_coeffs(x0, x1, x2)
    return c0 * f0 + c1 * f1 + c2 * f2
end

"""
    central_difference(f0, f1, f2, x0, x1, x2)
Given three points x0, x1, and x2, and the function values at those points, f0, f1, f2,
compute the second-order approximation of the derivative at x1

```jldoctest;setup = :(using HallThruster: central_difference)
f(x) = x^4
x0, x1, x2 = 1.9999999, 2, 2.0000001
central_difference(f(x0), f(x1), f(x2), x0, x1, x2)

# output

32.0
```
"""
@inline function central_difference(f0, f1, f2, x0, x1, x2)
    c0, c1, c2 = central_diff_coeffs(x0, x1, x2)
    return c0 * f0 + c1 * f1 + c2 * f2
end

"""
    backward_difference(f0, f1, f2, x0, x1, x2)
Given three points x0, x1, and x2, and the function values at those points, f0, f1, f2,
compute the second-order approximation of the derivative at x2

```jldoctest;setup = :(using HallThruster: backward_difference)
f(x) = x^4
x0, x1, x2 = 1.9999998, 1.9999999, 2
backward_difference(f(x0), f(x1), f(x2), x0, x1, x2)

# output

32.00000002980232
```
"""
@inline function backward_difference(f0, f1, f2, x0, x1, x2)
    c0, c1, c2 = backward_diff_coeffs(x0, x1, x2)
    return c0 * f0 + c1 * f1 + c2 * f2
end

"""
    second_deriv_central_diff(f0, f1, f2, x0, x1, x2)
Given three points x0, x1, and x2, and the function values at those points, f0, f1, f2,
compute the second-order approximation of the second derivative at x1

```jldoctest;setup = :(using HallThruster: second_deriv_central_diff)
f(x) = x^4
x0, x1, x2 = 1.999999, 2.0, 2.000001
second_deriv_central_diff(f(x0), f(x1), f(x2), x0, x1, x2)

# output

48.0
```
"""
@inline function second_deriv_central_diff(f0, f1, f2, x0, x1, x2)
    c0, c1, c2 = second_deriv_coeffs(x0, x1, x2)
    return c0 * f0 + c1 * f1 + c2 * f2
end