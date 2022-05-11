abstract type TransitionFunction end

Base.@kwdef struct SmoothIf <: TransitionFunction
    transition_length::Float64
end

(s::SmoothIf)(x, cutoff, y1, y2) = ((y2 - y1)*tanh((x-cutoff)/s.transition_length) + y1 + y2) / 2

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
