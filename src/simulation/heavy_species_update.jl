function iterate_heavy_species!(dU, U, params, reconstruct, source_heavy_species)
    (; cache, grid, ion_wall_losses, fluid_containers) = params

    # Compute edges and apply convective update
    _from_state_vector!(fluid_containers, U)

    update_convective_terms!(fluid_containers, grid, reconstruct, cache.dlnA_dz)
    source_heavy_species(fluid_containers, params)
    apply_reactions!(params.fluid_array, params)

    apply_ion_acceleration!(fluid_containers.isothermal, grid, cache)

    if ion_wall_losses
        apply_ion_wall_losses!(fluid_containers, params)
    end

    # Update maximum allowable timestep
    CFL = params.simulation.CFL
    min_dt_u = params.fluid_containers.continuity[1].max_timestep[]
    for fluid in params.fluid_containers.isothermal
        min_dt_u = min(min_dt_u, fluid.max_timestep[])
    end

    cache.dt[] = min(
        CFL * cache.dt_iz[],
        sqrt(CFL) * cache.dt_E[],
        CFL * min_dt_u,
    )

    _to_state_vector!(U, fluid_containers)

    # Update d/dt
    @. @views dU[1, :] = fluid_containers.continuity[1].dens_ddt

    for (i, fluid) in enumerate(fluid_containers.isothermal)
        @. @views dU[2 * i, :] = fluid.dens_ddt
        @. @views dU[2 * i + 1, :] = fluid.mom_ddt
    end

    return
end

"""
$(SIGNATURES)

Step the heavy species forward in time using the Strong-Stability-preserving RK22 (SSPRK22) algorithm.
This method is better known as Heun's method (https://en.wikipedia.org/wiki/Heun%27s_method).
The Butcher tableau of this method is

 0 │
 1 │ 1
 ─-╀─────────
   │ 1/2  1/2

The canonical form goes as follows for an ODE dy/dt = f(t, y) and step size h:

k_{n1} = f(t, y_n)
y_{n1} = y_n + h * k_{n1}
k_{n2} = f(t + h, y_{n1})
y_{n+1} = y_n + 0.5 * h * (k_{n1} + k_{n2})

As written, this requries three intermediate storage variables: k_{n1}, k_{n2}, and y_{n1}
We can reduce this to two using the following rearrangement

y_{n+1} = y_n / 2 + (y_n + h * k_{n1}) / 2 + h * k_{n2}
        = (y_n + y_{n1} + h * k_{n2}) / 2

With this, we do not need to store k_{n1} and k_{n2} separately and can instead reuse the same memory.
"""
function integrate_heavy_species!(u, params, reconstruct, source, dt)
    (; k, u1) = params.cache
    # First step
    # Compute slope k_{n1}
    iterate_heavy_species!(k, u, params, reconstruct, source)

    # Second step
    # Compute slope k_{n2}, resuing memory of k_{n1}
    @. u1 = u + dt * k
    stage_limiter!(u1, params)
    iterate_heavy_species!(k, u1, params, reconstruct, source)

    # Final step
    @. u = 0.5 * (u + u1 + dt * k)
    stage_limiter!(u, params)
    return
end

function update_heavy_species!(params)
    (; cache, propellants, propellants) = params

    # TODO: multiple propellants
    prop = propellants[1]
    mdot_a = prop.flow_rate_kg_s
    Ti = prop.ion_temperature_K

    # Apply fluid boundary conditions
    apply_left_boundary!(params.fluid_containers, cache, Ti, mdot_a, params.ingestion_density, params.anode_bc)
    apply_right_boundary!(params.fluid_containers)

    # Update ion variables as seen by electrons
    update_heavy_species_cache!(params.fluid_containers, params.cache, params.landmark)

    return
end

function update_heavy_species_cache!(fluids, cache, landmark)
    (; nn, ne, ni, ui, niui, Z_eff, ji, ϵ, nϵ, K) = cache

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

    # Compute electron mean energy for reactions
    @. ϵ = nϵ / ne
    if !landmark
        @. ϵ += K
    end

    return
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
        fluid.density[end] = fluid.density[end - 1]
    end

    @inbounds for fluid in fluids.isothermal
        fluid.density[end] = fluid.density[end - 1]
        fluid.momentum[end] = fluid.momentum[end - 1]
    end

    return
end
