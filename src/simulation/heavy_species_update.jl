function integrate_heavy_species!(params, user_source, dt)
    (;continuity, isothermal) = params.fluid_containers

    fluids = params.fluid_containers

    # We have temporarily replaced SSPRK22 with forward euler using two substeps.
    # This has the same cost and is about as stable
    num_substeps = 2
    _dt = dt / num_substeps
    for _ in 1:num_substeps
        # Update RHS of heavy species
        iterate_heavy_species!(fluids, params, user_source)

        # Integrate forward in time
        for fluid in continuity
            @. fluid.density += _dt * fluid.dens_ddt
        end

        for fluid in isothermal
            @. fluid.density += _dt * fluid.dens_ddt
            @. fluid.momentum += _dt * fluid.mom_ddt
        end

        # Apply stage limiter and check for infs or nans
        any_nan = stage_limiter!(continuity, isothermal)

        if any_nan
            return true
        end
    end

    # Apply fluid boundary conditions
    Ti = params.ion_temperature_K
    mdot_a = params.anode_mass_flow_rate
    ingestion_density = params.ingestion_density
    anode_bc = params.anode_bc

    apply_left_boundary!(params.fluid_containers, params.cache, Ti, mdot_a, ingestion_density, anode_bc)
    apply_right_boundary!(params.fluid_containers)

    # Update arrays in cache and apply boundaries
    update_heavy_species_cache!(params.fluid_containers, params.cache)

    return false
end

function iterate_heavy_species!(fluids, params, user_source!)
    (; cache, grid, ion_wall_losses) = params
    (; continuity, isothermal) = fluids

    # Compute edge fluxes and apply convective update
    update_convective_terms!(continuity, isothermal, grid, params.reconstruct, cache.dlnA_dz)

    apply_ion_acceleration!(isothermal, grid, cache)

    user_source!(fluids, params)

    if ion_wall_losses
        apply_ion_wall_losses!(continuity, isothermal, params)
    end

    apply_reactions!(params.fluid_arr, params)

    # Update maximum allowable timestep
    CFL = params.simulation.CFL
    min_dt_u = fluids.continuity[1].max_timestep[]
    for fluid in fluids.isothermal
        min_dt_u = min(min_dt_u, fluid.max_timestep[])
    end

    cache.dt[] = min(
        CFL * cache.dt_iz[],
        sqrt(CFL) * cache.dt_E[],
        CFL * min_dt_u,
    )

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
        reconstruct,
        dlnA_dz
    )

    for fluid in continuity
        compute_edge_states_continuity!(fluid, reconstruct)
        compute_fluxes_continuity!(fluid, grid)
        update_convective_terms_continuity!(fluid, grid)
    end

    for fluid in isothermal
        compute_edge_states_isothermal!(fluid, reconstruct)
        compute_fluxes_isothermal!(fluid, grid)
        update_convective_terms_isothermal!(fluid, grid, dlnA_dz)
    end

    return
end

function update_heavy_species_cache!(fluids, cache)
    (; nn, ne, ni, ui, niui, Z_eff, ji) = cache

    # Compute neutral number density
    @inbounds for fluid in fluids.continuity
        inv_m = inv(fluid.species.element.m)
        @. nn = fluid.density * inv_m
    end

    @. ne = 0
    @. Z_eff = 0
    @. ji = 0

    # Update plasma quantities
    @inbounds for (f, fluid) in enumerate(fluids.isothermal)
        inv_m = inv(fluid.species.element.m)
        Z = fluid.species.Z

        for i in eachindex(fluid.density)
            _ni = fluid.density[i] * inv_m
            _niui = fluid.momentum[i] * inv_m
            ni[f, i] = _ni
            niui[f, i] = _niui
            ui[f, i] = _niui / _ni
            ne[i] += Z * _ni
            Z_eff[i] += _ni
            ji[i] += Z * e * _niui
        end
    end

    @inbounds for i in eachindex(Z_eff)
        # Effective ion charge state (density-weighted average charge state)
        Z_eff[i] = max(1.0, ne[i] / Z_eff[i])
    end

    return
end

function stage_limiter!(continuity, isothermal)
    @inbounds for fluid in continuity
        if any(!isfinite, fluid.density)
            return true
        end

        min_density = MIN_NUMBER_DENSITY * fluid.species.element.m
        @simd for i in eachindex(fluid.density)
            fluid.density[i] = max(fluid.density[i], min_density)
        end
    end

    @inbounds for fluid in isothermal
        if any(!isfinite, fluid.density) || any(!isfinite, fluid.momentum)
            return true
        end

        m = fluid.species.element.m
        min_density = MIN_NUMBER_DENSITY * m
        @simd for i in eachindex(fluid.density)
            dens = fluid.density[i]
            vel = fluid.momentum[i] / dens
            fluid.density[i] = max(dens, min_density)
            fluid.momentum[i] = fluid.density[i] * vel
        end
    end
    return false
end

function apply_left_boundary!(fluids, cache, Ti, mdot_a, ingestion_density, anode_bc)
    Te_L = cache.Tev[1]

    bohm_factor = if anode_bc == :sheath
        Vs = cache.Vs[]
        # Compute sheath potential
        electron_repelling_sheath = Vs > 0
        if electron_repelling_sheath
            # Ion attracting/electron-repelling sheath, ions in pre-sheath attain reduced Bohm speed
            Vs_norm = Vs / Te_L
            # Compute correction factor (see Hara, PSST 28 (2019))
            χ = exp(-Vs_norm) / √(π * Vs_norm) / (1 + myerf(sqrt(Vs_norm)))
            inv(√(1 + χ))
        else
            # Ion-repelling sheath, ions have zero velocity at anode
            0.0
        end
    else
        1.0
    end

    # Add inlet neutral density
    # Add ingested mass flow rate at anode
    un = fluids.continuity[1].const_velocity
    neutral_density = mdot_a / cache.channel_area[1] / un
    neutral_density += ingestion_density

    @inbounds for fluid in fluids.isothermal

        mi = fluid.species.element.m
        Z = fluid.species.Z

        interior_density = fluid.density[2]
        interior_flux = fluid.momentum[2]
        interior_velocity = interior_flux / interior_density

        sound_speed = sqrt((kB * Ti + Z * e * Te_L) / mi)  # Ion acoustic speed
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

        neutral_density -= boundary_flux / un
        fluid.density[1] = boundary_density
        fluid.momentum[1] = boundary_flux
    end

    fluids.continuity[1].density[1] = neutral_density

    return
end

function apply_right_boundary!(fluids)
    @inbounds for fluid in fluids.continuity
        fluid.density[end] = fluid.density[end-1]
    end

    @inbounds for fluid in fluids.isothermal
        fluid.density[end] = fluid.density[end-1]
        fluid.momentum[end] = fluid.momentum[end-1]
    end

    return
end