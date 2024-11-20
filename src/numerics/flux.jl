"""
    $(TYPEDEF)
Defines a numerical flux function used in the heavy-species solve.
See `HallThruster.flux_functions` to see a list of available flux functions.
"""
struct FluxFunction{F}
    flux::F
end

(f::FluxFunction)(args...; kwargs...) = f.flux(args...; kwargs...)

# dummy methods for each flux function
function __rusanov end
function __global_lax_friedrichs end
function __HLLE end

# Define all flux functions
const rusanov = FluxFunction(__rusanov)
const global_lax_friedrichs = FluxFunction(__global_lax_friedrichs)
const HLLE = FluxFunction(__HLLE)

#=============================================================================
 Serialization
==============================================================================#
const flux_functions = (; rusanov, global_lax_friedrichs, HLLE)
Serialization.SType(::Type{T}) where {T <: FluxFunction} = Serialization.Enum()
Serialization.options(::Type{T}) where {T <: FluxFunction} = flux_functions

# Compute flux vector for 1D flows
@inline function flux(U::NTuple{1, T}, fluid) where {T}
    ρ = U[1]
    u = fluid.u
    return (ρ * u,)
end

# Compute flux vector for 2D flows
@inline function flux(U::NTuple{2, T}, fluid) where {T}
    _, ρu = U
    p = pressure(U, fluid)
    return (ρu, U[2]^2 / U[1] + p)
end

# Compute flux vector for 3D
@inline function flux(U::NTuple{3, T}, fluid) where {T}
    ρ, ρu, ρE = U
    u = ρu / ρ
    p = pressure(U, fluid)
    ρH = ρE + p
    return (ρu, U[2]^2 / U[1] + p, ρH * u)
end

# use fun metaprogramming create specialized flux versions for each type of fluid
for NUM_CONSERVATIVE in 1:3
    eval(
        quote
        @inbounds @fastmath function __rusanov(
                UL::NTuple{$NUM_CONSERVATIVE, T},
                UR::NTuple{$NUM_CONSERVATIVE, T},
                fluid,
                args...,
        ) where {T}
            γ = fluid.species.element.γ
            Z = fluid.species.Z

            uL = velocity(UL, fluid)
            uR = velocity(UR, fluid)
            TL = temperature(UL, fluid)
            TR = temperature(UR, fluid)
            aL = sound_speed(UL, fluid)
            aR = sound_speed(UL, fluid)

            sL_max = max(abs(uL - aL), abs(uL + aL))
            sR_max = max(abs(uR - aR), abs(uR + aR))

            smax = max(sL_max, sR_max)

            FL = flux(UL, fluid)
            FR = flux(UR, fluid)

            return @NTuple [0.5 * ((FL[j] + FR[j]) - smax * (UR[j] - UL[j]))
                            for
                            j in 1:($(NUM_CONSERVATIVE))]
        end

        @fastmath function __HLLE(
                UL::NTuple{$NUM_CONSERVATIVE, T},
                UR::NTuple{$NUM_CONSERVATIVE, T},
                fluid,
                args...,
        ) where {T}
            γ = fluid.species.element.γ
            Z = fluid.species.Z

            uL = velocity(UL, fluid)
            uR = velocity(UR, fluid)
            TL = temperature(UL, fluid)
            TR = temperature(UR, fluid)
            aL = sound_speed(UL, fluid)
            aR = sound_speed(UL, fluid)

            sL_min, sL_max = min(0, uL - aL), max(0, uL + aL)
            sR_min, sR_max = min(0, uR - aR), max(0, uR + aR)

            smin = min(sL_min, sR_min)
            smax = max(sL_max, sR_max)

            FL = flux(UL, fluid)
            FR = flux(UR, fluid)

            return @NTuple[0.5 * (FL[j] + FR[j]) -
                           0.5 * (smax + smin) / (smax - smin) * (FR[j] - FL[j]) +
                           smax * smin / (smax - smin) * (UR[j] - UL[j])
                           for
                           j in 1:($(NUM_CONSERVATIVE))]
        end

        @fastmath function __global_lax_friedrichs(
                UL::NTuple{$NUM_CONSERVATIVE, T},
                UR::NTuple{$NUM_CONSERVATIVE, T},
                fluid,
                max_wave_speed = 0.0,
                args...,
        ) where {T}
            γ = fluid.species.element.γ
            Z = fluid.species.Z

            uL = velocity(UL, fluid)
            uR = velocity(UR, fluid)

            FL = flux(UL, fluid)
            FR = flux(UR, fluid)

            return @NTuple[0.5 * (FL[j] + FR[j]) + 0.5 * max_wave_speed * (UL[j] - UR[j])
                           for
                           j in 1:($(NUM_CONSERVATIVE))]
        end
    end,
    )
end
