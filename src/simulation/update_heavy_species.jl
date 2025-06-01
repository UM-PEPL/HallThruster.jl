function iterate_heavy_species!(dU, U, params, scheme, user_source!)
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
        fluid.mom_ddt[i] = (fluid.flux_mom[left] - fluid.flux_mom[right]) / Δz - ρiui^2 / ρi * dlnA_dz[i]
    end

    return
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

function integrate_heavy_species!(U, params, scheme::HyperbolicScheme, user_source, dt)
    (; k) = params.cache

    # First step of SSPRK22
    iterate_heavy_species!(k, U, params, scheme, user_source)
    @. U = U + dt * k
    stage_limiter!(U, params)

    # Update arrays in cache
    update_heavy_species!(U, params)
    return
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

    @. ϵ = nϵ / ne + landmark * K
    return
end

function update_heavy_species!(U, params)
    (; index, grid, cache, mi, ncharge, landmark) = params

    # Apply fluid boundary conditions
    @views left_boundary_state!(U[:, 1], U, params)
    @views right_boundary_state!(U[:, end], U, params)

    update_heavy_species!(
        U, cache, index, grid.cell_centers, ncharge, mi, landmark,
    )
    return
end

function stage_limiter!(U, params)
    (; grid, index, min_Te, cache, mi, ncharge) = params
    stage_limiter!(U, grid.cell_centers, cache.nϵ, index, min_Te, ncharge, mi)
    return
end

function stage_limiter!(U, z_cell, nϵ, index, min_Te, ncharge, mi)
    min_density = MIN_NUMBER_DENSITY * mi
    @inbounds for i in eachindex(z_cell)
        U[index.ρn, i] = max(U[index.ρn, i], min_density)

        for Z in 1:ncharge
            density_floor = max(U[index.ρi[Z], i], min_density)
            velocity = U[index.ρiui[Z], i] / U[index.ρi[Z], i]
            U[index.ρi[Z], i] = density_floor
            U[index.ρiui[Z], i] = density_floor * velocity
        end
        nϵ[i] = max(nϵ[i], 1.5 * MIN_NUMBER_DENSITY * min_Te)
    end
    return
end

function left_boundary_state!(bc_state, U, params)
    index = params.index
    ncharge = params.ncharge
    mi = params.mi
    Ti = params.ion_temperature_K
    un = params.neutral_velocity
    mdot_a = params.anode_mass_flow_rate
    anode_bc = params.anode_bc
    ingestion_density = params.ingestion_density

    return left_boundary_state!(
        bc_state, U, index, ncharge, params.cache, mi,
        Ti, un, ingestion_density, mdot_a, anode_bc,
    )
end

function left_boundary_state!(
        bc_state, U, index, ncharge, cache, mi, Ti, un,
        ingestion_density, mdot_a, anode_bc,
    )
    if anode_bc == :sheath
        Vs = cache.Vs[]
        # Compute sheath potential
        electron_repelling_sheath = Vs > 0
        if electron_repelling_sheath
            # Ion attracting/electron-repelling sheath, ions in pre-sheath attain reduced Bohm speed
            Vs_norm = Vs / cache.Tev[1]
            # Compute correction factor (see Hara, PSST 28 (2019))
            χ = exp(-Vs_norm) / √(π * Vs_norm) / (1 + myerf(sqrt(Vs_norm)))
            bohm_factor = inv(√(1 + χ))
        else
            # Ion-repelling sheath, ions have zero velocity at anode
            bohm_factor = 0.0
        end
    else
        bohm_factor = 1.0
    end

    # Add inlet neutral density
    bc_state[index.ρn] = mdot_a / cache.channel_area[1] / un

    # Add ingested mass flow rate at anode
    bc_state[index.ρn] += ingestion_density

    return @inbounds for Z in 1:ncharge
        interior_density = U[index.ρi[Z], 2]
        interior_flux = U[index.ρiui[Z], 2]
        interior_velocity = interior_flux / interior_density

        sound_speed = sqrt((kB * Ti + Z * e * cache.Tev[1]) / mi)  # Ion acoustic speed
        boundary_velocity = -bohm_factor * sound_speed # Want to drive flow to (negative) bohm velocity

        if interior_velocity <= -sound_speed
            # Supersonic outflow, pure Neumann boundary condition
            boundary_density = interior_density
            boundary_flux = interior_flux
        else
            # Subsonic outflow, need to drive the flow toward sonic
            # For the isothermal Euler equations, the Riemann invariants are
            # J⁺ = u + c ln ρ
            # J⁻ = u - c ln ρ
            # For the boundary condition, we take c = u_bohm

            # 1. Compute outgoing characteristic using interior state
            J⁻ = interior_velocity - sound_speed * log(interior_density)

            # 2. Compute incoming characteristic using J⁻ invariant
            J⁺ = 2 * boundary_velocity - J⁻

            # 3. Compute boundary density using J⁺ and J⁻ invariants
            boundary_density = exp(0.5 * (J⁺ - J⁻) / sound_speed)

            # Compute boundary flux
            boundary_flux = boundary_velocity * boundary_density
        end

        bc_state[index.ρn] -= boundary_flux / un
        bc_state[index.ρi[Z]] = boundary_density
        bc_state[index.ρiui[Z]] = boundary_flux
    end
end

function right_boundary_state!(bc_state, U, params)
    (; index, ncharge) = params
    # Use Neumann boundary conditions for all neutral fluids
    bc_state[index.ρn] = U[index.ρn, end - 1]

    @inbounds for Z in 1:ncharge
        boundary_density = U[index.ρi[Z], end - 1]
        bc_state[index.ρi[Z]] = boundary_density        # Neumann BC for ion density at right boundary
        bc_state[index.ρiui[Z]] = U[index.ρiui[Z], end - 1] # Neumann BC for ion flux at right boundary
    end
    return
end
