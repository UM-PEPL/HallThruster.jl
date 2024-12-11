abstract type CurrentController end

struct NoController <: CurrentController end

@kwdef mutable struct PIDController <: CurrentController
    target_value::Float64
    proportional_constant::Float64 = 1.0
    integral_constant::Float64 = 0.0
    derivative_constant::Float64 = 0.0
    errors::NTuple{3, Float64} = (0.0, 0.0, 0.0)
    smoothing_frequency::Float64 = 0.0
    smoothed_value::Float64 = 0.0
end

#=============================================================================
 Serialization
==============================================================================#
Serialization.SType(::Type{T}) where {T <: CurrentController} = Serialization.TaggedUnion()
function Serialization.options(::Type{T}) where {T <: CurrentController}
    (; NoController, PIDController)
end
Serialization.exclude(::Type{T}) where {T <: CurrentController} = (:errors, :smoothed_val)

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

    Take derivatives again and rearrange: 

        u(t) = u(t-dt) +
            K_p (e(t) - e(t-dt)) + 
            K_i e(t) dt +
            K_d (e(t) - 2 e(t-dt) _ e(t-2dt))/dt.
    =#
    K_smooth = exp(-dt * pid.smoothing_frequency)
    pid.smoothed_value = K_smooth * present_value + (1 - K_smooth) * pid.smoothed_value

    pid.errors = (pid.target_value - present_value, pid.errors[1], pid.errors[2])

    P = pid.proportional_constant * (pid.errors[1] - pid.errors[2])
    I = pid.integral_constant * pid.errors[1] * dt
    D = pid.derivative_constant * (pid.errors[1] - 2 * pid.errors[2] + pid.errors[3]) / dt

    return control_value + P + I + D
end
