function iterate_heavy_species!(dU, U, params, scheme, sources; apply_boundary_conditions)
    (; index, cache, grid, simulation, ncharge, ion_wall_losses) = params
    (; source_neutrals, source_ion_continuity, source_ion_momentum) = sources
    (; F, UL, UR) = cache

    # Compute edges and apply convective update
    compute_fluxes!(F, UL, UR, U, params, scheme; apply_boundary_conditions)
    update_convective_term!(dU, U, F, grid, index, cache, ncharge)

    apply_user_ion_source_terms!(
        dU, U, params, source_neutrals, source_ion_continuity, source_ion_momentum,
    )

    apply_reactions!(dU, U, params)

    apply_ion_acceleration!(dU, U, params)

    if ion_wall_losses
        apply_ion_wall_losses!(dU, U, params)
    end

    update_timestep!(cache, dU, simulation.CFL, length(grid.cell_centers))

    return nothing
end

function update_timestep!(cache, dU, CFL, ncells)
    # Compute maximum allowable timestep
    cache.dt[] = min(
        CFL * cache.dt_iz[],                          # max ionization timestep
        sqrt(CFL) * cache.dt_E[],                     # max acceleration timestep
        CFL * minimum(@views cache.dt_u[1:(ncells - 1)]), # max fluid timestep
    )

    @. @views dU[:, 1] = 0.0
    return @. @views dU[:, end] = 0.0
end

function update_convective_term!(dU, U, F, grid, index, cache, ncharge)
    (; dA_dz, channel_area) = cache
    ncells = length(grid.cell_centers)
    return @inbounds for i in 2:(ncells - 1)
        left = left_edge(i)
        right = right_edge(i)

        Δz = grid.dz_cell[i]

        dlnA_dz = dA_dz[i] / channel_area[i]

        # Neutral flux
        dU[index.ρn, i] = (F[index.ρn, left] - F[index.ρn, right]) / Δz

        # Handle ions
        for Z in 1:ncharge
            # Ion fluxes
            # ∂ρ/∂t + ∂/∂z(ρu) = Q - ρu * ∂/∂z(lnA)
            ρi = U[index.ρi[Z], i]
            ρiui = U[index.ρiui[Z], i]
            dU[index.ρi[Z], i] = (F[index.ρi[Z], left] - F[index.ρi[Z], right]) / Δz - ρiui * dlnA_dz
            dU[index.ρiui[Z], i] = (F[index.ρiui[Z], left] - F[index.ρiui[Z], right]) / Δz - ρiui^2 / ρi * dlnA_dz
        end
    end
end

# Perform one step of the Strong-stability-preserving RK22 algorithm with the ion fluid
function integrate_heavy_species!(
        U, params, config::Config, dt, apply_boundary_conditions = true,
    )
    (; source_neutrals, source_ion_continuity, source_ion_momentum, scheme) = config
    sources = (; source_neutrals, source_ion_continuity, source_ion_momentum)

    return integrate_heavy_species!(U, params, scheme, sources, dt, apply_boundary_conditions)
end

function integrate_heavy_species!(
        U, params, scheme::HyperbolicScheme, sources, dt, apply_boundary_conditions = true,
    )
    (; k, u1) = params.cache
    # First step of SSPRK22
    iterate_heavy_species!(k, U, params, scheme, sources; apply_boundary_conditions)
    @. u1 = U + dt * k
    stage_limiter!(u1, params)

    # Second step of SSPRK22
    iterate_heavy_species!(k, u1, params, scheme, sources; apply_boundary_conditions)
    @. U = (U + u1 + dt * k) / 2
    return stage_limiter!(U, params)
end

function update_heavy_species!(U, cache, index, z_cell, ncharge, mi, landmark)
    (; nn, ne, ni, ui, niui, Z_eff, ji, K, ϵ, nϵ) = cache
    # Compute neutral number density
    inv_m = inv(mi)
    @. @views nn = U[index.ρn, :] * inv_m

    # Update plasma quantities
    @inbounds for i in eachindex(z_cell)
        # Compute ion derived quantities
        ne[i] = 0.0
        Z_eff[i] = 0.0
        ji[i] = 0.0
        @inbounds for Z in 1:(ncharge)
            _ni = U[index.ρi[Z], i] * inv_m
            _niui = U[index.ρiui[Z], i] * inv_m
            ni[Z, i] = _ni
            niui[Z, i] = _niui
            ui[Z, i] = _niui / _ni
            ne[i] += Z * _ni
            Z_eff[i] += _ni
            ji[i] += Z * e * _niui
        end

        # Effective ion charge state (density-weighted average charge state)
        Z_eff[i] = max(1.0, ne[i] / Z_eff[i])
    end

    return @. ϵ = nϵ / ne + landmark * K
end

function update_heavy_species!(U, params)
    (; index, grid, cache, mi, ncharge, landmark) = params

    # Apply fluid boundary conditions
    @views left_boundary_state!(U[:, 1], U, params)
    @views right_boundary_state!(U[:, end], U, params)

    return update_heavy_species!(
        U, cache, index, grid.cell_centers, ncharge, mi, landmark,
    )
end
