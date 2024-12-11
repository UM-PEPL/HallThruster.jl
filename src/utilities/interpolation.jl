struct LinearInterpolation{X, Y}
    xs::X
    ys::Y
    function LinearInterpolation(x, y; resample_uniform = false, resample_factor = 2)
        if length(x) != length(y)
            throw(ArgumentError("x and y must have same length"))
        end

        # Resample x-values to be uniform to speed up finding indices
        if resample_uniform
            xmin, xmax = extrema(x)
            num = length(x) * resample_factor
            resampled_x = range(xmin, xmax, num)
            resampled_y = [interpolate(_x, x, y) for _x in resampled_x]
            return LinearInterpolation(resampled_x, resampled_y; resample_uniform = false)
        else
            return new{typeof(x), typeof(y)}(x, y)
        end
    end
end

function find_left_index(value, array)
    N = length(array)

    left = 0
    right = N + 1

    @inbounds while (right - left) > 1
        mid = (left + right) >>> 0x01

        cond = array[mid] > value
        not_cond = !cond

        # conditional assignments
        right = cond * mid + not_cond * right
        left = not_cond * mid + cond * left
    end
    return left
end

@fastmath function interpolate(x::T, xs, ys; use_log = false) where {T}
    i = find_left_index(x, xs)
    i < 1 && return ys[1] / oneunit(T)
    i ≥ length(xs) && return ys[end] / oneunit(T)
    use_log && return exp(lerp(x, xs[i], xs[i + 1], log(ys[i]), log(ys[i + 1])))
    return lerp(x, xs[i], xs[i + 1], ys[i], ys[i + 1])
end

@inbounds @fastmath function interpolate(x::T, xs::StepRangeLen, ys; use_log) where {T}
    dx_inv = inv(xs.step.hi)
    i = 1 + floor(Int, (x - xs[1]) * dx_inv)
    i < 1 && return T(ys[1])
    i ≥ length(xs) && return T(ys[end])
    t = (x - xs[i]) * dx_inv
    if (use_log)
        y1 = log(ys[i + 1])
        y0 = log(ys[i])
        return exp(muladd(t, y1 - y0, y0))
    else
        return muladd(t, ys[i + 1] - ys[i], ys[i])
    end
end

function (itp::LinearInterpolation)(x::T; use_log = false) where {T}
    interpolate(x, itp.xs, itp.ys; use_log)
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
