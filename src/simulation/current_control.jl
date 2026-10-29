abstract type CurrentController end

struct NoController <: CurrentController end

@kwdef mutable struct PIDController <: CurrentController
    target_value::Float64
    proportional_constant::Float64 = 1.0
    integral_constant::Float64 = 0.0
    derivative_constant::Float64 = 0.0
    errors::NTuple{3, Float64} = (0.0, 0.0, 0.0)
    previous_dt::Float64 = 0.0
    previous_error_slope::Float64 = 0.0
    smoothing_frequency::Float64 = 0.0
    smoothed_value::Float64 = 0.0
end

#=============================================================================
 Serialization
==============================================================================#
Serialization.SType(::Type{T}) where {T <: CurrentController} = Serialization.TaggedUnion()
function Serialization.options(::Type{T}) where {T <: CurrentController}
    return (; NoController, PIDController)
end
Serialization.exclude(::Type{T}) where {T <: CurrentController} =
    (:errors, :previous_dt, :previous_error_slope, :smoothed_value)

#=============================================================================
 Definitions
==============================================================================#

apply_controller(::CurrentController, ::Any, control_value, ::Any) = control_value

function apply_controller(pid::PIDController, present_value, control_value, dt)
    #=
    Use discrete form of PID control
    derived from
    https://en.wikipedia.org/wiki/Proportional%E2%80%93integral%E2%80%93derivative_controller#Discrete_implementation

    # Nomenclature:
    - e(t): target_value - present_value(t)
    - u(t): control_value
    - K_p:  proportional constant
    - K_i = K_p / T_i:  integral constant
    - K_d = K_p T_d:    derivative constant

    # Derivation

    Start with basic PID equation:

        u(t) = K_p e(t) + K_i ∫ e(t) dt + K_d e'(t).

    Take the derivative of both sides with respect to time:

         u'(t) = K_p e'(t) + K_i e(t) + K_d e''(t).

    Discretize derivatives (f'(t) = (f(t) - f(t-dt))/dt):

        u(t) - u(t-dt) = K_p (e(t) - e(t-dt)) + K_i e(t) dt + K_d (e'(t) - e'(t-dt)).

    For unequal timesteps h and h_prev, the derivative contribution is

        u(t) = u(t-dt) +
            K_p (e(t) - e(t-dt)) + 
            K_i e(t) dt +
            K_d ((e(t) - e(t-h))/h - (e(t-h) - e(t-h-h_prev))/h_prev).
    =#
    isfinite(dt) && dt > 0 || throw(ArgumentError("Controller timestep must be positive and finite."))
    pid.smoothed_value = if pid.smoothing_frequency > 0
        K_smooth = exp(-dt * pid.smoothing_frequency)
        K_smooth * pid.smoothed_value + (-expm1(-dt * pid.smoothing_frequency)) * present_value
    else
        present_value
    end

    pid.errors = (pid.target_value - pid.smoothed_value, pid.errors[1], pid.errors[2])

    P = pid.proportional_constant * (pid.errors[1] - pid.errors[2])
    I = pid.integral_constant * pid.errors[1] * dt
    # The first measurement has no preceding sample. Seed the derivative at
    # zero instead of differentiating the fictitious initial zero error over
    # the tiny setup interval. Later updates use their actual sample spacing.
    present_slope = pid.previous_dt > 0 ?
        (pid.errors[1] - pid.errors[2]) / dt : 0.0
    D = pid.derivative_constant * (present_slope - pid.previous_error_slope)
    pid.previous_error_slope = present_slope
    pid.previous_dt = dt

    return control_value + P + I + D
end
