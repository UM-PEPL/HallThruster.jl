@inline function reconstruct(uⱼ₋₁, uⱼ, uⱼ₊₁, limiter)
    r = (uⱼ₊₁ - uⱼ) / (uⱼ - uⱼ₋₁)
    Δu = 0.25 * limiter(r) * (uⱼ₊₁ - uⱼ₋₁)
    return uⱼ - Δu, uⱼ + Δu
end

function compute_edge_states!(UL, UR, U, params, config; apply_boundary_conditions = false)
    (; is_velocity_index, index, ncharge) = params
    (; scheme) = config

    compute_edge_states!(UL, UR, U, scheme.limiter, is_velocity_index, scheme.reconstruct)

    if apply_boundary_conditions
        @views left_boundary_state!(UL[:, 1], U, params, config)
        @views right_boundary_state!(UR[:, end], U, index, ncharge)
    else
        @. @views UL[:, 1] = U[:, 1]
        @. @views UR[:, end] = U[:, end]
    end
end

function compute_edge_states!(UL, UR, U, limiter, is_velocity_index, do_reconstruct)
    (nvars, ncells) = size(U)

    # compute left and right edge states
    if do_reconstruct
        @inbounds for j in 1:nvars
            if is_velocity_index[j] # reconstruct velocity as primitive variable instead of momentum density
                for i in 2:(ncells - 1)
                    u₋ = U[j, i - 1] / U[j - 1, i - 1]
                    uᵢ = U[j, i] / U[j - 1, i]
                    u₊ = U[j, i + 1] / U[j - 1, i + 1]
                    uR, uL = reconstruct(u₋, uᵢ, u₊, limiter)

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

                    UR[j, left_edge(i)], UL[j, right_edge(i)] = reconstruct(
                        u₋, uᵢ, u₊, limiter,)
                end
            end
        end
    else
        @inbounds for i in 2:(ncells - 1), j in 1:nvars
            UL[j, right_edge(i)] = U[j, i]
            UR[j, left_edge(i)]  = U[j, i]
        end
    end
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
end

function compute_fluxes!(F, UL, UR, flux_function, λ_global, grid, fluids, index, ncharge)
    @inbounds for i in eachindex(grid.edges)
        # Neutral fluxes at edge i
        left_state_n  = (UL[index.ρn, i],)
        right_state_n = (UR[index.ρn, i],)

        F[index.ρn, i] = flux_function(
            left_state_n, right_state_n, fluids[1], λ_global[1],)[1]

        # Ion fluxes at edge i
        for Z in 1:ncharge
            left_state_i        = (UL[index.ρi[Z], i], UL[index.ρiui[Z], i])
            right_state_i       = (UR[index.ρi[Z], i], UR[index.ρiui[Z], i])
            fluid_ind           = Z + 1
            F_mass, F_momentum  = flux_function(left_state_i, right_state_i, fluids[fluid_ind], λ_global[fluid_ind])
            F[index.ρi[Z], i]   = F_mass
            F[index.ρiui[Z], i] = F_momentum
        end
    end
end

function compute_fluxes!(F, UL, UR, U, params, config; apply_boundary_conditions = false)
    (; index, fluids, grid, cache) = params
    (; λ_global, dt_u) = cache
    (; scheme, ncharge) = config

    # Reconstruct the states at the left and right edges using MUSCL scheme
    compute_edge_states!(UL, UR, U, params, config; apply_boundary_conditions)
    compute_wave_speeds!(λ_global, dt_u, UL, UR, U, grid, fluids, index, ncharge)
    compute_fluxes!(F, UL, UR, scheme.flux_function, λ_global, grid, fluids, index, ncharge)

    return F
end
