@inline check_r(r) = isfinite(r) && r >= 0
@inline van_leer_limiter(r) = check_r(r) * (4r / (r + 1)^2)
@inline primitive_velocity(momentum, density) = density > 0 ? momentum / density : 0.0

"""Refresh one fluid's primitive-velocity cache from its conservative state."""
function update_primitive_velocity!(fluid)
    @inbounds @simd for i in eachindex(fluid.vel_prim)
        fluid.vel_prim[i] = primitive_velocity(fluid.momentum[i], fluid.density[i])
    end
    return nothing
end

@inline function reconstruct(uⱼ₋₁, uⱼ, uⱼ₊₁)
    Δu_L = uⱼ - uⱼ₋₁
    Δu_R = uⱼ₊₁ - uⱼ
    same_sign = signbit(Δu_L) == signbit(Δu_R)
    valid = same_sign && Δu_L != 0 && isfinite(Δu_L) && isfinite(Δu_R)
    # Harmonic form of the Van Leer slope: r / (1 + r)^2 * (Δu_L + Δu_R).
    Δu = valid ? Δu_L * Δu_R / (Δu_L + Δu_R) : 0.0
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

    # Get edge density from boundary cells
    # More details in `compute_edge_states_isothermal`
    fluid.dens_L[1] = fluid.density[1]
    fluid.dens_R[1] = fluid.density[2]
    fluid.dens_L[end] = fluid.density[end - 1]
    fluid.dens_R[end] = fluid.density[end]
    return
end

function compute_edge_states_isothermal!(fluid, do_reconstruct)
    (; density, vel_prim, dens_L, dens_R, vel_L, vel_R) = fluid
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
            u₋ = vel_prim[i - 1]
            uᵢ = vel_prim[i]
            u₊ = vel_prim[i + 1]
            uR, uL = reconstruct(u₋, uᵢ, u₊)
            vel_L[iR] = uL
            vel_R[iL] = uR
        end
    else
        @inbounds for i in 2:(N - 1)
            iL, iR = left_edge(i), right_edge(i)
            dens_L[iR] = density[i]
            dens_R[iL] = density[i]
            velocity = vel_prim[i]
            vel_L[iR] = velocity
            vel_R[iL] = velocity
        end
    end

    #        Left boundary
    #             V   First
    #            =|   cell    |
    #       o----=|-----o-----|---
    # Cell: 1    =|     2     |
    # Edge:     L 1 R       L 2 R
    #
    # We calculate boundary properties at edge 1,
    # and extrapolate linearly to set the cell 1 properties.
    # The flux at the left edge should be set according to the edge property,
    # which we can compute by averaging the cell 1 and cell 2 properties.

    fluid.dens_L[1] = fluid.density[1]
    fluid.dens_R[1] = fluid.density[2]
    fluid.dens_L[end] = fluid.density[end - 1]
    fluid.dens_R[end] = fluid.density[end]

    fluid.vel_L[1] = vel_prim[1]
    fluid.vel_R[1] = vel_prim[2]
    fluid.vel_L[end] = vel_prim[end - 1]
    fluid.vel_R[end] = vel_prim[end]

    return
end

function compute_fluxes_continuity!(fluid, grid)
    (; flux_dens, dens_L, dens_R, wave_speed, const_velocity) = fluid
    smax = wave_speed[]

    # The neutral wave speed and grid never change during a simulation, so its
    # CFL limit only needs to be computed on the first flux update.
    if fluid.max_timestep[] <= 0
        min_timestep = Inf
        @inbounds for i in eachindex(grid.dz_edge)
            min_timestep = min(min_timestep, grid.dz_edge[i] / smax)
        end
        fluid.max_timestep[] = min_timestep
    end

    return @inbounds for i in eachindex(fluid.dens_L)
        ρ_L, ρ_R = dens_L[i], dens_R[i]
        flux_dens[i] = 0.5 * (const_velocity * (ρ_L + ρ_R) - smax * (ρ_R - ρ_L))
    end
end

function compute_fluxes_isothermal!(fluid, grid)
    (; flux_dens, flux_mom, dens_L, dens_R, vel_L, vel_R) = fluid
    a = fluid.sound_speed
    RT = a^2 / fluid.species.element.γ

    min_timestep = Inf

    @inbounds for i in eachindex(dens_L)
        ρ_L, ρ_R = dens_L[i], dens_R[i]
        u_L, u_R = vel_L[i], vel_R[i]
        ρu_L, ρu_R = ρ_L * u_L, ρ_R * u_R

        # For nonnegative sound speed, max(|u - a|, |u + a|) = |u| + a.
        smax = max(abs(u_L), abs(u_R)) + a
        min_timestep = min(min_timestep, grid.dz_edge[i] / smax)

        flux_mom_L = ρ_L * (u_L^2 + RT)
        flux_mom_R = ρ_R * (u_R^2 + RT)

        flux_dens[i] = 0.5 * ((ρu_L + ρu_R) - smax * (ρ_R - ρ_L))
        flux_mom[i] = 0.5 * ((flux_mom_L + flux_mom_R) - smax * (ρu_R - ρu_L))
    end

    fluid.max_timestep[] = min_timestep
    return
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
        fluid.mom_ddt[i] = (fluid.flux_mom[left] - fluid.flux_mom[right]) / Δz -
            ρiui * fluid.vel_prim[i] * dlnA_dz[i]
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
        update_primitive_velocity!(fluid)
        compute_edge_states_isothermal!(fluid, reconstruct)
        compute_fluxes_isothermal!(fluid, grid)
        update_convective_terms_isothermal!(fluid, grid, dlnA_dz)
    end

    return
end
