@inline function linear_transition(x, cutoff, L, y1, y2)
    x1 = cutoff - L / 2
    x2 = cutoff + L / 2
    if x < x1
        return y1
    elseif x > x2
        return y2
    else
        return lerp(x, x1, x2, y1, y2)
    end
end

function smooth_if(x, cutoff, y1, y2, L)
    return ((y2 - y1) * tanh((x - cutoff) / (L / 4)) + y1 + y2) / 2
end

function smooth!(x, x_cache; iters = 1)
    if iters > 0
        x_cache .= x
        x_cache[1] = x[2]
        x_cache[end - 1] = x[end]
        for i in 2:(length(x) - 1)
            if i == 2 || i == length(x) - 1
                x_cache[i] = 0.5 * x[i] + 0.25 * (x[i - 1] + x[i + 1])
            else
                x_cache[i] = 0.4 * x[i] + 0.24 * (x[i - 1] + x[i + 1]) +
                             0.06 * (x[i - 2] + x[i + 2])
            end
        end
        x .= x_cache
        smooth!(x, x_cache; iters = iters - 1)
    else
        return x
    end
end
