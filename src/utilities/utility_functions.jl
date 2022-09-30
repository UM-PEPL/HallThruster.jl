struct LinearInterpolation{X<:Number,Y<:Number}
    xs::Vector{X}
    ys::Vector{Y}
    function LinearInterpolation(x, y)
        if length(x) != length(y)
            throw(ArgumentError("x and y must have same length"))
        else
            return new{typeof(x[1]),typeof(y[1])}(x, y)
        end
    end
end

function interpolate(x::T, xs, ys; use_log = false) where {T}
    if x ≤ xs[1]
        return ys[1] / oneunit(T)
    elseif x ≥ xs[end]
        return ys[end] / oneunit(T)
    end
    i = find_left_index(x, xs)
    itp = if use_log
        exp(lerp(x, xs[i], xs[i+1], log(ys[i]), log(ys[i+1])))
    else
        lerp(x, xs[i], xs[i+1], ys[i], ys[i+1])
    end
    return itp
end

function (itp::LinearInterpolation)(x::T; use_log = false) where {T}
    xs, ys = itp.xs, itp.ys
    interpolate(x, xs, ys; use_log)
end

"""
    lerp(x, x0, x1, y0, y1)
Interpolate between two points (x0, y0) and (x1, y1)
```jldoctest;setup = :(using HallThruster: lerp)
julia> lerp(0.5, 0.0, 1.0, 0.0, 2.0)
1.0
"""
@inline function lerp(x, x0, x1, y0, y1)
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
end


function find_left_index(value, array)
    N = length(array)

    #=if value ≥ array[end]
        return N
    elseif value < array[1]
        return 0
    elseif value == array[1]
        return 1
    end=#

    left = 0
    right = N

    @inbounds while (right - left) > 1
        mid = (left + right) >>> 0x01
        if array[mid] > value
            right = mid
        else
            left = mid
        end
    end
    return left
end


left_edge(i) = i - 1
right_edge(i) = i

@inline electron_density(U, p, i) = sum(Z * U[p.index.ρi[Z], i] for Z in 1:p.config.ncharge) / p.config.propellant.m

@inline inlet_neutral_density(config) = config.anode_mass_flow_rate / config.neutral_velocity / config.thruster.geometry.channel_area

@inline background_neutral_density(config) = config.propellant.m * config.background_pressure / kB / config.background_neutral_temperature

@inline background_neutral_velocity(config) = -sqrt(kB * config.background_neutral_temperature / config.propellant.m)

@inline ion_current_density(U, p, i) = sum(Z * e * U[p.index.ρiui[Z], i] for Z in 1:p.config.ncharge) / p.config.propellant.m

function discharge_current(U::Matrix, params)
    (;A_ch, cache, Δz_edge, ϕ_L, ϕ_R) = params
    (;∇pe, μ, ne, ji, Vs) = cache

    ncells = size(U, 2)

    int1 = 0.0
    int2 = 0.0

    @inbounds for i in 1:ncells-1
        Δz = Δz_edge[i]

        int1_1 = (ji[i] / e / μ[i] + ∇pe[i]) / ne[i]
        int1_2 = (ji[i+1] / e / μ[i+1] + ∇pe[i+1]) / ne[i+1]

        int1 += 0.5 * Δz * (int1_1 + int1_2)

        int2_1 = inv(e * ne[i] * μ[i])
        int2_2 = inv(e * ne[i+1] * μ[i+1])

        int2 += 0.5 * Δz * (int2_1 + int2_2)
    end

    ΔV = ϕ_L + Vs[] - ϕ_R

    I = A_ch * (ΔV + int1) / int2

    return I
end

function cumtrapz(x, y, y0 = zero(typeof(y[1] * x[1])))
    int = zeros(typeof(y[1] * x[1]), length(x))
    cumtrapz!(int, x, y, y0)
end

function cumtrapz!(cache, x, y, y0 = zero(typeof(y[1] * x[1])))

    cache[1] = y0
    @inbounds for i in 2:lastindex(x)
        Δx = x[i] - x[i-1]
        cache[i] = cache[i-1] + 0.5 * Δx * (y[i] + y[i-1])
    end

    return cache
end

function tridiagonal_forward_sweep!(A::Tridiagonal, b)
    n = length(A.d)

    @inbounds for i in 2:n
        w = A.dl[i - 1] / A.d[i - 1]
        A.d[i] = A.d[i] - w * A.du[i - 1]
        b[i] = b[i] - w * b[i - 1]
    end
end

function tridiagonal_backward_sweep!(y, A::Tridiagonal, b)
    n = length(A.d)
    y[n] = b[n] / A.d[n]
    # back-substitution
    @inbounds for i in (n - 1):-1:1
        y[i] = (b[i] - A.du[i] * y[i + 1]) / A.d[i]
    end
end

# our matrix is diagonally dominant so we can use Thomas' algorithm to solve
# the tridiagonal system
function tridiagonal_solve!(y, A, b)
    tridiagonal_forward_sweep!(A, b)
    return tridiagonal_backward_sweep!(y, A, b)
end

function tridiagonal_solve(A, b)
    y = similar(b)
    A′ = copy(A)
    b′ = copy(b)
    tridiagonal_solve!(y, A′, b′)
    return y
end

function myerf(x)
    x < 0 && return -myerf(-x)
    x_sqrt_pi = x * √(π)
    x_squared = x^2

    h = x_sqrt_pi + (π-2)*x_squared
    g = h / (1 + h)
    return 1 - exp(-x_squared) / x_sqrt_pi * g
end
