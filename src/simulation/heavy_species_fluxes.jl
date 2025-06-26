function stage_limiter!(U, params)
    (; grid, index, min_Te, cache, propellants) = params

    # TODO: multiple propellants + fluid containers
    ncharge = propellants[1].max_charge
    mi = propellants[1].gas.m
    return stage_limiter!(U, grid.cell_centers, cache.nϵ, index, min_Te, ncharge, mi)
end

function stage_limiter!(U, z_cell, nϵ, index, min_Te, ncharge, mi)
    min_density = MIN_NUMBER_DENSITY * mi
    return @inbounds for i in eachindex(z_cell)
        U[index.ρn, i] = max(U[index.ρn, i], min_density)

        for Z in 1:ncharge
            density_floor = max(U[index.ρi[Z], i], min_density)
            velocity = U[index.ρiui[Z], i] / U[index.ρi[Z], i]
            U[index.ρi[Z], i] = density_floor
            U[index.ρiui[Z], i] = density_floor * velocity
        end
        nϵ[i] = max(nϵ[i], 1.5 * MIN_NUMBER_DENSITY * min_Te)
    end
end

"""
    $(TYPEDEF)
Defines a numerical flux function used in the heavy-species solve.
See `HallThruster.flux_functions` to see a list of available flux functions.
"""
struct FluxFunction{F}
    flux::F
end

(f::FluxFunction)(args...; kwargs...) = f.flux(args...; kwargs...)

# declare flux functions
function __rusanov end

# Define all flux functions
const rusanov = FluxFunction(__rusanov)

#=============================================================================
 Serialization
==============================================================================#
const flux_functions = (; rusanov)
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
                uL = velocity(UL, fluid)
                uR = velocity(UR, fluid)
                aL = sound_speed(UL, fluid)
                aR = sound_speed(UL, fluid)

                sL_max = max(abs(uL - aL), abs(uL + aL), abs(uL))
                sR_max = max(abs(uR - aR), abs(uR + aR), abs(uR))

                smax = max(sL_max, sR_max)

                FL = flux(UL, fluid)
                FR = flux(UR, fluid)

                return @NTuple [
                    0.5 * ((FL[j] + FR[j]) - smax * (UR[j] - UL[j]))
                        for
                        j in 1:($(NUM_CONSERVATIVE))
                ]
            end
        end,
    )
end

@inline check_r(r) = isfinite(r) && r >= 0
@inline van_leer_limiter(r) = check_r(r) * (4r / (r + 1)^2)

@inline function reconstruct(uⱼ₋₁, uⱼ, uⱼ₊₁)
    r = (uⱼ₊₁ - uⱼ) / (uⱼ - uⱼ₋₁)
    Δu = 0.25 * van_leer_limiter(r) * (uⱼ₊₁ - uⱼ₋₁)
    return uⱼ - Δu, uⱼ + Δu
end

function compute_edge_states!(UL, UR, U, params, do_reconstruct; apply_boundary_conditions = false)
    (nvars, ncells) = size(U)

    # compute left and right edge states
    if do_reconstruct
        @inbounds for j in 1:nvars
            if params.is_velocity_index[j] # reconstruct velocity as primitive variable instead of momentum density
                for i in 2:(ncells - 1)
                    u₋ = U[j, i - 1] / U[j - 1, i - 1]
                    uᵢ = U[j, i] / U[j - 1, i]
                    u₊ = U[j, i + 1] / U[j - 1, i + 1]
                    uR, uL = reconstruct(u₋, uᵢ, u₊)

                    ρL = UL[j - 1, right_edge(i)] #use previously-reconstructed edge density to compute momentum
                    ρR = UR[j - 1, left_edge(i)]
                    UL[j, right_edge(i)] = uL * ρL
                    UR[j, left_edge(i)] = uR * ρR
                end
            else
                for i in 2:(ncells - 1)
                    u₋ = U[j, i - 1]
                    uᵢ = U[j, i]
                    u₊ = U[j, i + 1]

                    UR[j, left_edge(i)], UL[j, right_edge(i)] = reconstruct(u₋, uᵢ, u₊)
                end
            end
        end
    else
        @inbounds for i in 2:(ncells - 1), j in 1:nvars
            UL[j, right_edge(i)] = U[j, i]
            UR[j, left_edge(i)] = U[j, i]
        end
    end

    if apply_boundary_conditions
        @views left_boundary_state!(UL[:, 1], U, params)
        @views right_boundary_state!(UR[:, end], U, params)
    else
        @. @views UL[:, 1] = U[:, 1]
        @. @views UR[:, end] = U[:, end]
    end
    return
end

function compute_wave_speeds!(λ_global, dt_u, UL, UR, U, grid, fluids, index, ncharge)
    # Compute maximum wave speed in domain and use this to update the max allowable timestep, if using adaptive timestepping
    @inbounds for i in eachindex(grid.edges)
        # Compute wave speeds for each component of the state vector.
        # The only wave speed for neutrals is the neutral convection velocity
        neutral_fluid = fluids[1]
        U_neutrals = (U[index.ρn, i],)
        u = velocity(U_neutrals, neutral_fluid)
        λ_global[1] = abs(u)

        # Ion wave speeds
        for Z in 1:ncharge
            fluid_ind = Z + 1
            fluid = fluids[fluid_ind]
            γ = fluid.species.element.γ
            UL_ions = (UL[index.ρi[Z], i], UL[index.ρiui[Z], i])
            UR_ions = (UR[index.ρi[Z], i], UR[index.ρiui[Z], i])

            uL = velocity(UL_ions, fluid)
            uR = velocity(UR_ions, fluid)
            aL = sound_speed(UL_ions, fluid)
            aR = sound_speed(UL_ions, fluid)

            # Maximum wave speed
            s_max = max(abs(uL + aL), abs(uL - aL), abs(uR + aR), abs(uR - aR))

            # a Δt / Δx = 1 for CFL condition, user-supplied CFL number restriction applied later, in update_values
            dt_max = grid.dz_edge[i] / s_max
            dt_u[i] = dt_max

            # Update maximum wavespeeds and maximum allowable timestep
            λ_global[fluid_ind] = max(s_max, λ_global[fluid_ind])
        end
    end
    return
end

function compute_fluxes!(F, UL, UR, λ_global, grid, fluids, index, ncharge)
    @inbounds for i in eachindex(grid.edges)
        # Neutral fluxes at edge i
        left_state_n = (UL[index.ρn, i],)
        right_state_n = (UR[index.ρn, i],)

        F[index.ρn, i] = rusanov(left_state_n, right_state_n, fluids[1], λ_global[1])[1]

        # Ion fluxes at edge i
        for Z in 1:ncharge
            left_state_i = (UL[index.ρi[Z], i], UL[index.ρiui[Z], i])
            right_state_i = (UR[index.ρi[Z], i], UR[index.ρiui[Z], i])
            fluid_ind = Z + 1
            F_mass, F_momentum = rusanov(left_state_i, right_state_i, fluids[fluid_ind], λ_global[fluid_ind])
            F[index.ρi[Z], i] = F_mass
            F[index.ρiui[Z], i] = F_momentum
        end
    end
    return
end

function compute_fluxes!(F, UL, UR, U, params, reconstruct; apply_boundary_conditions = false)
    (; index, fluids, grid, cache, propellants) = params
    (; λ_global, dt_u) = cache

    ncharge = propellants[1].max_charge

    # Reconstruct the states at the left and right edges using MUSCL scheme
    compute_edge_states!(UL, UR, U, params, reconstruct; apply_boundary_conditions)
    compute_wave_speeds!(λ_global, dt_u, UL, UR, U, grid, fluids, index, ncharge)
    compute_fluxes!(F, UL, UR, λ_global, grid, fluids, index, ncharge)

    return F
end
