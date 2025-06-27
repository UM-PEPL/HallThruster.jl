function integrate_heavy_species!(fluid_containers, params, user_source, dt)
    # Do one timestep forward, returning `true` if we found a NaN or Inf
    step_heavy_species!(fluid_containers, params, user_source, dt) && return true
    # Update properties that interface with electrons
    update_heavy_species!(params)
    return false
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
function step_heavy_species!(fluid_containers, params, source, dt)
    # First step
    # Compute slope k_{n1}
    compute_heavy_species_derivatives!(fluid_containers, params, source)

    # Copy density and momentum to dens_cache and mom_cache for all fluids
    # and update density and momentum to y_{n1}
    for fluid in fluid_containers.continuity
        @. fluid.dens_cache = fluid.density
        @. fluid.density += dt * fluid.dens_ddt
    end

    for fluid in fluid_containers.isothermal
        @. fluid.dens_cache = fluid.density
        @. fluid.mom_cache = fluid.momentum
        @. fluid.density += dt * fluid.dens_ddt
        @. fluid.momentum += dt * fluid.mom_ddt
    end

    # Apply stage limiter, returning true if a NaN or an Inf is detected
    stage_limiter!(fluid_containers) && return true

    # Second step
    # Compute slope k_{n2}, reusing memory of k_{n1}
    compute_heavy_species_derivatives!(fluid_containers, params, source)

    # Final step
    @inbounds for fluid in fluid_containers.continuity
        @. fluid.density = 0.5 * (fluid.density + fluid.dens_cache + dt * fluid.dens_ddt)
    end

    @inbounds for fluid in fluid_containers.isothermal
        @. fluid.density = 0.5 * (fluid.density + fluid.dens_cache + dt * fluid.dens_ddt)
        @. fluid.momentum = 0.5 * (fluid.momentum + fluid.mom_cache + dt * fluid.mom_ddt)
    end

    # Apply stage limiter, returning true if a NaN or an Inf is detected
    stage_limiter!(fluid_containers) && return true

    return false
end

# Populate dens_ddt and mom_ddt for all fluid containers
function compute_heavy_species_derivatives!(fluid_containers, params, source_heavy_species)
    (; cache, grid, ion_wall_losses, reconstruct) = params

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

    return
end

function stage_limiter!(fluid_containers)
    @inbounds for fluid in fluid_containers.continuity
        if any(!isfinite, fluid.density)
            return true
        end

        min_density = MIN_NUMBER_DENSITY * fluid.species.element.m
        @simd for i in eachindex(fluid.density)
            fluid.density[i] = max(fluid.density[i], min_density)
        end
    end

    @inbounds for fluid in fluid_containers.isothermal
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

function update_heavy_species!(params)
    (; cache, propellants, anode_bc, ingestion_density, total_flow_rate_kg_s) = params

    # Apply left boundary conditions per-propellant
    for (propellant, fluids) in zip(propellants, params.fluids_by_propellant)
        apply_left_boundary!(fluids, propellant, cache, anode_bc, ingestion_density, total_flow_rate_kg_s)
    end

    # Apply right boundary conditions for all propellants
    apply_right_boundary!(params.fluid_containers)

    # Update ion variables as seen by electrons
    update_heavy_species_cache!(params.fluid_containers, cache, params.landmark)

    return
end

function update_heavy_species_cache!(fluids, cache, landmark)
    (; nn, ne, ni, ui, niui, Z_eff, ji, ϵ, nϵ, K, m_eff) = cache

    # Compute neutral number density
    @inbounds for fluid in fluids.continuity
        inv_m = inv(fluid.species.element.m)
        @. nn = fluid.density * inv_m
    end

    @. ne = 0
    @. Z_eff = 0
    @. ji = 0
    @. m_eff = 0

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
            ji[i] += Z * e * _niui
            # First pass, store total ion number density in Z_eff
            Z_eff[i] += _ni
            m_eff[i] += fluid.density[i]
        end
    end

    @inbounds for i in eachindex(Z_eff)
        # Effective ion charge state (density-weighted average charge state)
        inv_ion_number_density = inv(Z_eff[i])
        Z_eff[i] = ne[i] * inv_ion_number_density
        m_eff[i] = m_eff[i] * inv_ion_number_density
    end

    # Compute electron mean energy for reactions
    @. ϵ = nϵ / ne
    if !landmark
        @. ϵ += K
    end

    return
end

function apply_left_boundary!(fluids, propellant, cache, anode_bc, ingestion_density, total_flow_rate_kg_s)
    Te_L = cache.Tev[1]
    Ti = propellant.ion_temperature_K
    mdot_a = propellant.flow_rate_kg_s

    # Scale ingested neutral density by ratio of the flow rate of this propellant to the total flow rate
    ingestion_density = ingestion_density * mdot_a / total_flow_rate_kg_s

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
    un = fluids.continuity[].const_velocity
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

    fluids.continuity[].density[1] = neutral_density

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
