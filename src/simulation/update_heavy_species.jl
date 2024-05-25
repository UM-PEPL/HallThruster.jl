function iterate_heavy_species!(dU, U, params; apply_boundary_conditions = true)
    (;index, Δz_cell, config, cache, ncells) = params
    (;
        source_neutrals, source_ion_continuity, source_ion_momentum,
        ncharge,
    ) = config

    (;F, UL, UR, channel_area, dA_dz) = cache

    compute_fluxes!(F, UL, UR, U, params; apply_boundary_conditions)
    @inbounds for i in 2:ncells-1
        left = left_edge(i)
        right = right_edge(i)

        Δz = Δz_cell[i]

        dlnA_dz = dA_dz[i] / channel_area[i]

        # Neutral flux
        dU[index.ρn, i] = (F[index.ρn, left] - F[index.ρn, right]) / Δz

        # User-provided neutral source term
        dU[index.ρn, i] += source_neutrals[1](U, params, i)

        # Handle ions
        for Z in 1:ncharge
            # Ion fluxes
            # ∂ρ/∂t + ∂/∂z(ρu) = Q - ρu * ∂/∂z(lnA)
            ρi = U[index.ρi[Z], i]
            ρiui = U[index.ρiui[Z], i]
            dU[index.ρi[Z]  , i] = (F[index.ρi[Z],   left] - F[index.ρi[Z],   right]) / Δz - ρiui * dlnA_dz
            dU[index.ρiui[Z], i] = (F[index.ρiui[Z], left] - F[index.ρiui[Z], right]) / Δz - ρiui^2 / ρi * dlnA_dz

            # User-provided ion source terms
            dU[index.ρi[Z],   i] += source_ion_continuity[Z](U, params, i)
            dU[index.ρiui[Z], i] += source_ion_momentum[Z  ](U, params, i)
        end
    end

    @. @views dU[:, 1] = 0.0
    @. @views dU[:, end] = 0.0

    return nothing
end

function heavy_species_source_terms!(dU, U, params, dt)
    (;ncells, cache, config, CFL, index) = params
    (;ncharge) = config
    dU .= 0.0

    # update electron density in interior cells
    @inbounds for i in 1:ncells
        cache.ne[i] = 0.0
        for Z in 1:ncharge
            cache.ni[Z, i] = U[index.ρi[Z], i] / config.propellant.m
            cache.ne[i] += Z * cache.ni[Z, i]
        end
    end

    # reset timestep
    params.cache.dt[] = Inf

    @inbounds for i in 2:ncells-1
        apply_ion_acceleration!(dU, params, i)
        apply_reactions!(dU, U, params, i)

        if config.ion_wall_losses
            apply_ion_wall_losses!(dU, U, params, i)
        end

        # update allowable timestep
        cache.dt_cell[i] = min(
            sqrt(CFL) * cache.dt_E[i],
            CFL * cache.dt_iz[i],
            CFL * cache.dt_u[left_edge(i)],
            CFL * cache.dt_u[right_edge(i)]
        )

        cache.dt[] = min(cache.dt_cell[i], cache.dt[])
    end

    cache.dt_cell[1]   = cache.dt_cell[2]
    cache.dt_cell[end] = cache.dt_cell[end-1]

    return nothing
end

# Perform one step of the Strong-stability-preserving RK22 algorithm to the ion fluid
function integrate_heavy_species!(U, params, dt, apply_boundary_conditions = true)
    (;k, u1) = params.cache

    # First step of SSPRK22 for convective terms
    iterate_heavy_species!(k, U, params; apply_boundary_conditions)
    @. u1 = U + dt * k
    stage_limiter!(u1, params)

    # Second step of SSPRK22 for convective terms
    iterate_heavy_species!(k, u1, params; apply_boundary_conditions)
    @. U = (U + u1 + dt * k) / 2
    stage_limiter!(U, params)

    # Apply source terms
    heavy_species_source_terms!(k, U, params, params.dt[])
    @. U += k * dt

    return nothing
end

function update_heavy_species!(U, params)
    (;index, ncells, cache) = params
    (;nn, ne, ni, ui, niui, Z_eff, ji) = cache
    mi = params.config.propellant.m

    # Apply fluid boundary conditions
    @views left_boundary_state!(U[:, 1], U, params)
    @views right_boundary_state!(U[:, end], U, params)

    # Update plasma quantities
    @inbounds for i in 1:ncells
        # Compute number density for each neutral fluid
        nn[i] = U[index.ρn, i] / params.config.propellant.m

        # Compute ion derived quantities
        ne[i] = 0.0
        Z_eff[i] = 0.0
        ji[i] = 0.0
        @inbounds for Z in 1:params.config.ncharge
            _ni = U[index.ρi[Z], i] / mi
            _niui = U[index.ρiui[Z], i] / mi
            ni[Z, i] = _ni
            niui[Z, i] = _niui
            ui[Z, i] = _niui / _ni
            ne[i] += Z * _ni
            Z_eff[i] += _ni
            ji[i] += Z * e * _niui
        end

        # Compute electron number density, making sure it is above floor
        ne[i] = max(params.config.min_number_density, ne[i])

        # Effective ion charge state (density-weighted average charge state)
        Z_eff[i] = max(1.0, ne[i] / Z_eff[i])
    end
end
