@inline function reconstruct(uⱼ₋₁, uⱼ, uⱼ₊₁, limiter)
    r = (uⱼ₊₁ - uⱼ) / (uⱼ - uⱼ₋₁)
    Δu = 0.25 * limiter(r) * (uⱼ₊₁ - uⱼ₋₁)
    return uⱼ - Δu, uⱼ + Δu
end

function compute_edge_states!(fluid::ContinuityFluid, limiter, do_reconstruct)
    (; density, dens_L, dens_R) = fluid
    N = length(fluid.density)

    if do_reconstruct
        @inbounds for i in 2:(N - 1)
            iL, iR = left_edge(i), right_edge(i)

            # Reconstruct density
            u₋ = density[i - 1]
            uᵢ = density[i]
            u₊ = density[i + 1]
            dens_R[iL], dens_L[iR] = reconstruct(u₋, uᵢ, u₊, limiter)
        end
    else
        @inbounds for i in 2:(N - 1)
            iL, iR = left_edge(i), right_edge(i)
            dens_L[iR] = density[i]
            dens_R[iL] = density[i]
        end
    end

    fluid.dens_L[1] = fluid.density[1]
    return fluid.dens_R[end] = fluid.density[end]
end

function compute_edge_states!(fluid::IsothermalFluid, limiter, do_reconstruct)
    (; density, momentum, dens_L, dens_R, mom_L, mom_R) = fluid
    N = length(fluid.density)

    if do_reconstruct
        @inbounds for i in 2:(N - 1)
            iL, iR = left_edge(i), right_edge(i)

            # Reconstruct density
            u₋ = density[i - 1]
            uᵢ = density[i]
            u₊ = density[i + 1]
            dens_R[iL], dens_L[iR] = reconstruct(u₋, uᵢ, u₊, limiter)

            # Reconstruct velocity as primitive variable instead of momentum density
            u₋ = momentum[i - 1] / u₋
            uᵢ = momentum[i] / uᵢ
            u₊ = momentum[i + 1] / u₊
            uR, uL = reconstruct(u₋, uᵢ, u₊, limiter)
            mom_L[iR] = uL * dens_L[iR]
            mom_R[iL] = uR * dens_R[iL]
        end
    else
        @inbounds for i in 2:(N - 1)
            iL, iR = left_edge(i), right_edge(i)
            dens_L[iR] = density[i]
            dens_R[iL] = density[i]
            mom_L[iR] = momentum[i]
            mom_R[iL] = momentum[i]
        end
    end

    fluid.dens_L[1] = fluid.density[1]
    fluid.dens_R[end] = fluid.density[end]
    fluid.mom_L[1] = fluid.momentum[1]
    return fluid.mom_R[end] = fluid.momentum[end]
end

function compute_fluxes!(fluid::ContinuityFluid, grid)
    (; flux_dens, dens_L, dens_R, wave_speed, velocity) = fluid
    smax = wave_speed[]
    fluid.max_timestep[] = Inf

    @inbounds for i in eachindex(fluid.dens_L)
        ρ_L, ρ_R = dens_L[i], dens_R[i]
        flux_dens[i] = 0.5 * (velocity * (ρ_L + ρ_R) - smax * (ρ_R - ρ_L))
        fluid.max_timestep[] = min(fluid.max_timestep[], grid.dz_edge[i] / smax)
    end
end

function compute_fluxes!(fluid::IsothermalFluid, grid)
    (; flux_dens, flux_mom, dens_L, dens_R, mom_L, mom_R, wave_speed) = fluid
    a = fluid.sound_speed
    RT = a^2 / fluid.species.element.γ

    max_wave_speed = 0.0
    fluid.max_timestep[] = Inf

    @inbounds for i in eachindex(dens_L)
        ρ_L, ρ_R = dens_L[i], dens_R[i]
        ρu_L, ρu_R = mom_L[i], mom_R[i]

        u_L, u_R = ρu_L / ρ_L, ρu_R / ρ_R

        smax = max(abs(u_L - a), abs(u_L + a), abs(u_R - a), abs(u_R + a))
        fluid.max_timestep[] = min(fluid.max_timestep[], grid.dz_edge[i] / smax)
        max_wave_speed = max(smax, max_wave_speed)

        flux_mom_L = ρ_L * (u_L^2 + RT)
        flux_mom_R = ρ_R * (u_R^2 + RT)

        flux_dens[i] = 0.5 * ((ρu_L + ρu_R) - smax * (ρ_R - ρ_L))
        flux_mom[i] = 0.5 * ((flux_mom_L + flux_mom_R) - smax * (ρu_R - ρu_L))
    end

    return wave_speed[] = max_wave_speed
end