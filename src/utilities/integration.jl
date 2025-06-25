function cumtrapz(x, y, y0 = zero(typeof(y[1] * x[1])))
    int = zeros(typeof(y[1] * x[1]), length(x))
    return cumtrapz!(int, x, y, y0)
end

function cumtrapz!(cache, x, y, y0 = zero(typeof(y[1] * x[1])))
    cache[1] = y0
    @inbounds for i in 2:lastindex(x)
        Δx = x[i] - x[i - 1]
        cache[i] = cache[i - 1] + 0.5 * Δx * (y[i] + y[i - 1])
    end

    return cache
end
