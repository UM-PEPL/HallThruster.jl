abstract type CurrentController end

struct NoController <: CurrentController end

@kwdef mutable struct PIDController <: CurrentController
    target_value::Float64
    proportional_constant::Float64 = 1.0
    integral_constant::Float64 = 0.0
    derivative_constant::Float64 = 0.0
    smoothing_time::Float64 = 5e-4
    smoothed_value::Float64 = 0.0
    errors::NTuple{3, Float64} = (0.0, 0.0, 0.0)
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

function apply_pid_controller(pid::PIDController, present_value, control_value, dt)

    # Apply exponential smoothing to input signal
    K_smooth = 1 - exp(-dt / pid.smoothing_time)
    pid.smoothed_value = K_smooth * present_value + (1 - K_smooth) * pid.smoothed_value

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

    Take derivatives again:

        u(t) - u(t-dt) = 
            K_p (e(t) - e(t-dt)) + K_i e(t) dt + K_d ((e(t) - e(t-dt)) - (e(t-dt) - e(t-2dt)))/dt.

    Rearrange:

        u(t) = u(t-dt) + C1 e(t) + C2 e(t-dt) + C3 e(t-2dt),

    with
        C1 = K_p + K_i dt + K_d/dt,
        C2 = -K_p - 2 K_d/dt,
        C3 = K_d/dt.
    =#

    K_p = pid.proportional_constant
    K_i = pid.integral_constant
    K_d = pid.derivative_constant
    err = pid.target_value - pid.smoothed_value
    pid.errors = (err, pid.errors[1], pid.errors[2])

    C1 = K_p + K_i * dt + K_d / dt
    C2 = -K_p - 2 * K_d / dt
    C3 = K_d / dt

    return control_value + C1 * pid.errors[1] + C2 * pid.errors[2] + C3 * pid.errors[3]
end
