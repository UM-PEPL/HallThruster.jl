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

@inline check_r(r) = isfinite(r) && r >= 0
@inline van_leer_limiter(r) = check_r(r) * (4r / (r + 1)^2)

@inline function reconstruct(uⱼ₋₁, uⱼ, uⱼ₊₁)
    r = (uⱼ₊₁ - uⱼ) / (uⱼ - uⱼ₋₁)
    Δu = 0.25 * van_leer_limiter(r) * (uⱼ₊₁ - uⱼ₋₁)
    return uⱼ - Δu, uⱼ + Δu
end

function compute_edge_states_continuity!(fluid, do_reconstruct)
    (; density, dens_L, dens_R) = fluid
    N = length(fluid.density)

    if do_reconstruct
        @inbounds for i in 2:(N - 1)
            iL, iR = left_edge(i), right_edge(i)
            # Reconstruct density
            u₋ = density[i - 1]
            uᵢ = density[i]
            u₊ = density[i + 1]
            dens_R[iL], dens_L[iR] = reconstruct(u₋, uᵢ, u₊)
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

function compute_edge_states_isothermal!(fluid, do_reconstruct)
    (; density, momentum, dens_L, dens_R, mom_L, mom_R) = fluid
    N = length(fluid.density)

    if do_reconstruct
        @inbounds for i in 2:(N - 1)
            iL, iR = left_edge(i), right_edge(i)

            # Reconstruct density
            u₋ = density[i - 1]
            uᵢ = density[i]
            u₊ = density[i + 1]
            dens_R[iL], dens_L[iR] = reconstruct(u₋, uᵢ, u₊)

            # Reconstruct velocity as primitive variable instead of momentum density
            u₋ = momentum[i - 1] / u₋
            uᵢ = momentum[i] / uᵢ
            u₊ = momentum[i + 1] / u₊
            uR, uL = reconstruct(u₋, uᵢ, u₊)
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

function compute_fluxes_continuity!(fluid, grid)
    (; flux_dens, dens_L, dens_R, wave_speed, const_velocity) = fluid
    smax = wave_speed[]
    fluid.max_timestep[] = Inf

    return @inbounds for i in eachindex(fluid.dens_L)
        ρ_L, ρ_R = dens_L[i], dens_R[i]
        flux_dens[i] = 0.5 * (const_velocity * (ρ_L + ρ_R) - smax * (ρ_R - ρ_L))
        fluid.max_timestep[] = min(fluid.max_timestep[], grid.dz_edge[i] / smax)
    end
end

function compute_fluxes_isothermal!(fluid, grid)
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

function update_convective_terms_continuity!(fluid, grid)
    ncells = length(grid.cell_centers)
    @inbounds for i in 2:(ncells - 1)
        left, right = left_edge(i), right_edge(i)
        Δz = grid.dz_cell[i]
        fluid.dens_ddt[i] = (fluid.flux_dens[left] - fluid.flux_dens[right]) / Δz
    end
    return
end

function update_convective_terms_isothermal!(fluid, grid, dlnA_dz)
    ncells = length(grid.cell_centers)

    @inbounds for i in 2:(ncells - 1)
        left, right = left_edge(i), right_edge(i)
        Δz = grid.dz_cell[i]

        # ∂ρ/∂t + ∂/∂z(ρu) = Q - ρu * ∂/∂z(lnA)
        ρi = fluid.density[i]
        ρiui = fluid.momentum[i]
        fluid.dens_ddt[i] = (fluid.flux_dens[left] - fluid.flux_dens[right]) / Δz - ρiui * dlnA_dz[i]
        fluid.mom_ddt[i] = (fluid.flux_mom[left] - fluid.flux_mom[right]) / Δz - ρiui^2 / ρi * dlnA_dz[i]
    end

    return
end

function update_convective_terms!(fluid_containers, grid, reconstruct, dlnA_dz)

    for fluid in fluid_containers.continuity
        compute_edge_states_continuity!(fluid, reconstruct)
        compute_fluxes_continuity!(fluid, grid)
        update_convective_terms_continuity!(fluid, grid)
    end

    for fluid in fluid_containers.isothermal
        compute_edge_states_isothermal!(fluid, reconstruct)
        compute_fluxes_isothermal!(fluid, grid)
        update_convective_terms_isothermal!(fluid, grid, dlnA_dz)
    end

    return
end
