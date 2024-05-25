
function update_heavy_species!(dU, U, params, t; apply_boundary_conditions = true)

    (;index, Δz_cell, config, cache) = params
    (;
        source_neutrals, source_ion_continuity, source_ion_momentum,
        ncharge, ion_wall_losses
    ) = config

    (;F, UL, UR, channel_area, dA_dz) = cache

    ncells = size(U, 2)

    compute_fluxes!(F, UL, UR, U, params; apply_boundary_conditions)

    params.cache.dt[] = Inf

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

        apply_ion_acceleration!(dU, U, params, i)
        apply_reactions!(dU, U, params, i)

        if ion_wall_losses
            apply_ion_wall_losses!(dU, U, params, i)
        end

        params.cache.dt_cell[i] = min(
            sqrt(params.CFL) * params.cache.dt_E[i],
            params.CFL * params.cache.dt_iz[i],
            params.CFL * params.cache.dt_u[left],
            params.CFL * params.cache.dt_u[right]
        )

        params.cache.dt[] = min(params.cache.dt_cell[i], params.cache.dt[])
    end

    params.cache.dt_cell[1] = params.cache.dt_cell[2]
    params.cache.dt_cell[end] = params.cache.dt_cell[end-1]

    @. @views dU[:, 1] = 0.0
    @. @views dU[:, end] = 0.0

    return nothing
end
