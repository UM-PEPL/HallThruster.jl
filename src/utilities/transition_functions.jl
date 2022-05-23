abstract type TransitionFunction end

Base.@kwdef struct SmoothIf <: TransitionFunction
    transition_length::Float64
end

@inline smooth_if(x, cutoff, y1, y2, L) = ((y2 - y1)*tanh((x-cutoff)/(L/4)) + y1 + y2) / 2

(s::SmoothIf)(x, cutoff, y1, y2) = smooth_if(x, cutoff, y1, y2, s.transition_length)

Base.@kwdef struct SmoothStep <: TransitionFunction
    transition_length::Float64
end

@inline smoothstep(x, x1, x2, y1, y2) = let t = (x - x1) / (x2 - x1)
    y = if t < 0
        0.0
    elseif t > 1
        1.0
    else
        3 * t^2 - 2 * t^3
    end
    return (y2 - y1) * y + y1
end

(s::SmoothStep)(x, cutoff, y1, y2) = let L = s.transition_length
    smoothstep(x, cutoff - L/2, cutoff + L/2, y1, y2)
end

Base.@kwdef struct QuadraticTransition <: TransitionFunction
    transition_length::Float64
    offset::Float64 = 0.0
end

function (q::QuadraticTransition)(x, cutoff, y1, y2)
    x′ = x - q.offset
    if x′ < cutoff
        return y1
    else
        return y1 + (y2 - y1) * (x′ - cutoff)^2 / ((x - cutoff)^2 + q.transition_length^2)
    end
end

Base.@kwdef struct LinearTransition <: TransitionFunction
    transition_length::Float64
    offset::Float64 = 0.0
end

function(ℓ::LinearTransition)(x, cutoff, y1, y2)
    x′ = x - ℓ.offset
    L = ℓ.transition_length
    x1 = cutoff - L/2
    x2 = cutoff + L/2
    if x′ < x1
        return y1
    elseif x′ > x2
        return y2
    else
        return lerp(x, x1, x2, y1, y2)
    end
end

struct StepFunction <: TransitionFunction end

(::StepFunction)(x, cutoff, y1, y2) = x ≤ cutoff ? y1 : y2
