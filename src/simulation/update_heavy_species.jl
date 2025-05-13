function iterate_heavy_species!(dU, U, params, scheme, user_source!; apply_boundary_conditions)
    (; cache, grid, ion_wall_losses, fluid_containers) = params
    (; continuity, isothermal) = fluid_containers

    # Populate fluid containers to compute fluxes
    _from_state_vector!(continuity, isothermal, U)

    # Compute edge fluxes and apply convective update
    update_convective_terms!(continuity, isothermal, grid, scheme, cache.dlnA_dz)

    apply_ion_acceleration!(isothermal, grid, cache)

    user_source!(fluid_containers, params)

    if ion_wall_losses
        apply_ion_wall_losses!(continuity, isothermal, params)
    end

    apply_reactions!(params.fluid_arr, params)

    # Transfer fluid container d/dt to dU
    index = 1
    for fluid in continuity
        @. @views dU[index, :] = fluid.dens_ddt
        index += 1
    end

    for fluid in isothermal
        @. @views dU[index, :] = fluid.dens_ddt
        @. @views dU[index + 1, :] = fluid.mom_ddt
        index += 2
    end

    # Update maximum allowable timestep
    CFL = params.simulation.CFL
    min_dt_u = fluid_containers.continuity[1].max_timestep[]
    for fluid in fluid_containers.isothermal
        min_dt_u = min(min_dt_u, fluid.max_timestep[])
    end

    cache.dt[] = min(
        CFL * cache.dt_iz[],
        sqrt(CFL) * cache.dt_E[],
        CFL * min_dt_u,
    )

    # Set changes in left and right cells to zero
    @. @views dU[:, 1] = 0.0
    @. @views dU[:, end] = 0.0

    return nothing
end


function update_convective_terms_continuity!(fluid, grid)
    ncells = length(grid.cell_centers)
    return @inbounds for i in 2:(ncells - 1)
        left, right = left_edge(i), right_edge(i)
        Δz = grid.dz_cell[i]
        fluid.dens_ddt[i] = (fluid.flux_dens[left] - fluid.flux_dens[right]) / Δz
    end
end

function update_convective_terms_isothermal!(fluid, grid, dlnA_dz)
    ncells = length(grid.cell_centers)

    return @inbounds for i in 2:(ncells - 1)
        left, right = left_edge(i), right_edge(i)
        Δz = grid.dz_cell[i]

        # ∂ρ/∂t + ∂/∂z(ρu) = Q - ρu * ∂/∂z(lnA)
        ρi = fluid.density[i]
        ρiui = fluid.momentum[i]
        fluid.dens_ddt[i] = (fluid.flux_dens[left] - fluid.flux_dens[right]) / Δz - ρiui * dlnA_dz[i]
        fluid.mom_ddt[i] = (fluid.flux_mom[left] - fluid.flux_mom[right]) / Δz - ρiui^2 / ρi * dlnA_dz[i]
    end
end

function update_convective_terms!(
        continuity,
        isothermal,
        grid,
        scheme,
        dlnA_dz
    )

    for fluid in continuity
        compute_edge_states_continuity!(fluid, scheme.limiter, scheme.reconstruct)
        compute_fluxes_continuity!(fluid, grid)
        update_convective_terms_continuity!(fluid, grid)
    end

    for fluid in isothermal
        compute_edge_states_isothermal!(fluid, scheme.limiter, scheme.reconstruct)
        compute_fluxes_isothermal!(fluid, grid)
        update_convective_terms_isothermal!(fluid, grid, dlnA_dz)
    end


    return
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
        U, params, scheme::HyperbolicScheme, user_source, dt, apply_boundary_conditions = true,
    )
    (; k, u1) = params.cache

    # First step of SSPRK22
    iterate_heavy_species!(k, U, params, scheme, user_source; apply_boundary_conditions)
    @. u1 = U + dt * k
    stage_limiter!(u1, params)

    # Second step of SSPRK22
    iterate_heavy_species!(k, u1, params, scheme, user_source; apply_boundary_conditions)
    @. U = (U + u1 + dt * k) / 2
    stage_limiter!(U, params)

    # Update arrays in cache
    return update_heavy_species!(U, params)
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
