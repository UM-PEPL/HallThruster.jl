struct LinearInterpolation{X, Y}
    xs::X
    ys::Y
    function LinearInterpolation(x, y; resample_uniform = false, resample_factor = 2)
        if length(x) != length(y)
            throw(ArgumentError("x and y must have same length"))
        end

        # Resample xs to be a uniform LinRange to speed up finding indices
        if resample_uniform
            xmin, xmax = extrema(x)
            num = length(x) * resample_factor
            itp = LinearInterpolation(x, y)

            resampled_x = LinRangeWrapper(LinRange(xmin, xmax, num))
            resampled_y = itp.(resampled_x)
            return LinearInterpolation(resampled_x, resampled_y; resample_uniform = false)
        else
            return new{typeof(x),typeof(y)}(x, y)
        end
    end
end

struct LinRangeWrapper{T1, T2}
    r::LinRange{T1, T2}
    A::T1
    B::T1
end

Base.length(r::LinRangeWrapper) = length(r.r)
Base.getindex(r::LinRangeWrapper, args...) = getindex(r.r, args...)
Base.iterate(r::LinRangeWrapper, args...) = iterate(r.r, args...)
Base.firstindex(r::LinRangeWrapper, args...) = firstindex(r.r, args...)
Base.lastindex(r::LinRangeWrapper, args...) = lastindex(r.r, args...)

function LinRangeWrapper(r)
    A = r.len / (r.stop - r.start)
    b = -A * r.start
    return LinRangeWrapper(r, A, b)
end

@fastmath function interpolate(x::T, xs, ys; use_log = false) where {T}
    isnan(x) && return NaN
    i = find_left_index(x, xs)

    i < 1          && return ys[1] / oneunit(T)
    i ≥ length(xs) && return ys[end] / oneunit(T)

    use_log && return exp(lerp(x, xs[i], xs[i+1], log(ys[i]), log(ys[i+1])))
    return lerp(x, xs[i], xs[i+1], ys[i], ys[i+1])
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
    t = (x - x0) / (x1 - x0)
    return muladd(t, (y1 - y0), y0)
end

# ceil(Int, log2(x)) done with bitwise ops
function bitwise_log2ceil(x)
    count = 1
    while x != 1
        x >>= 0x01 # divide x by 2
        count +=1
    end
    return count
end

@inline function find_left_index(val, r::LinRangeWrapper)
    #ceil(Int, clamp(muladd(r.A, val, r.B), 0, r.r.len))
    return searchsortedfirst(r.r, val) - 1
end

function find_left_index(value, array)
    N = length(array)

    left = 0
    right = N+1

    @inbounds while (right - left) > 1
        mid = (left + right) >>> 0x01

        cond = array[mid] > value
        not_cond = !cond

        # conditional assignments
        right = cond * mid + not_cond * right
        left  = not_cond * mid + cond * left
    end
    return left
end

@inline left_edge(i) = i - 1
@inline right_edge(i) = i

@inline electron_density(U, p, i) = sum(Z * U[p.index.ρi[Z], i] for Z in 1:p.config.ncharge) / p.config.propellant.m

@inline inlet_neutral_density(config) = config.anode_mass_flow_rate / config.neutral_velocity / config.thruster.geometry.channel_area

@inline background_neutral_density(config) = config.propellant.m * config.background_pressure / kB / config.background_neutral_temperature

@inline background_neutral_velocity(config) = -sqrt(kB * config.background_neutral_temperature / config.propellant.m)

@inline ion_current_density(U, p, i) = sum(Z * e * U[p.index.ρiui[Z], i] for Z in 1:p.config.ncharge) / p.config.propellant.m

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
