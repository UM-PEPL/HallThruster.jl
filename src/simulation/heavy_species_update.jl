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
        for (_, fluids) in zip(propellants, params.fluids_by_propellant)
            apply_ion_wall_losses!(fluids, params)
        end
    end

    # Update maximum allowable timestep
    CFL = params.simulation.CFL
    min_dt_u = Inf
    for fluid in params.fluid_array
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

        min_density = MIN_NUMBER_DENSITY * fluid.species.element.m
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
    (; cache, propellants, anode_bc, ingestion_flow_rates) = params

    # Apply left boundary conditions per-propellant
    for (i, (propellant, fluids)) in enumerate(zip(propellants, params.fluids_by_propellant))
        apply_left_boundary!(fluids, propellant, cache, anode_bc, ingestion_flow_rates[i])
    end

    # Apply right boundary conditions for all propellants
    apply_right_boundary!(params.fluid_containers)

    # Update ion variables as seen by electrons
    update_heavy_species_cache!(params.fluid_containers, cache, params.landmark)

    return
end

function update_heavy_species_cache!(fluids, cache, landmark)
    (; nn, ne, Z_eff, ji, ϵ, nϵ, K, m_eff, avg_ion_vel, avg_neutral_vel) = cache

    @. ne = 0
    @. ji = 0
    @. m_eff = 0
    @. Z_eff = 0
    @. nn = 0
    @. avg_ion_vel = 0
    @. avg_neutral_vel = 0

    # Compute neutral number density
    # TODO: this computes total neutral number density, not per species
    @inbounds for fluid in fluids.continuity
        _nn = fluid.density / fluid.species.element.m
        @. nn += _nn
        @. avg_neutral_vel += _nn * fluid.const_velocity
    end


    # Update plasma quantities
    @inbounds for fluid in fluids.isothermal
        inv_m = inv(fluid.species.element.m)
        Z = fluid.species.Z

        for i in eachindex(fluid.density)
            _ni = fluid.density[i] * inv_m
            _niui = fluid.momentum[i] * inv_m
            ne[i] += Z * _ni
            ji[i] += Z * e * _niui
            avg_ion_vel[i] += _niui
            # First pass, store total ion mass density in m_eff and ion number density in Z_eff
            Z_eff[i] += _ni
            m_eff[i] += fluid.density[i]
        end
    end

    @. avg_neutral_vel /= nn

    @. ne = max(ne, MIN_NUMBER_DENSITY)

    # Inverse ion density
    @. Z_eff = inv(Z_eff)
    @. avg_ion_vel *= Z_eff
    @. m_eff *= Z_eff
    @. Z_eff = ne * Z_eff

    # Compute electron mean energy for reactions
    @. ϵ = nϵ / ne
    if !landmark
        @. ϵ += K
    end

    return
end

#===============================================================================
Boundary conditions
===============================================================================#

function apply_left_boundary!(fluids, propellant, cache, anode_bc, ingestion_flow_rate)
    Te_L = 0.5 * (cache.Tev[2] + cache.Tev[1])         # eV
    Ti = propellant.ion_temperature_K                 # K
    mdot_a = propellant.flow_rate_kg_s

    kTe_J = e * Te_L      # electron energy in Joules
    kTi_J = kB * Ti        # ion thermal energy in Joules
    γ = kTe_J / kTi_J

    # Sheath-edge electronegativity
    # Use first interior cell (index 2) as the sheath-edge estimate.
    n_pos_charge = 0.0
    n_neg_charge = 0.0
    @inbounds for fluid in fluids.isothermal
        Z = fluid.species.Z
        n = fluid.density[2]
        if Z > 0
            n_pos_charge += Z * n
        elseif Z < 0
            n_neg_charge += abs(Z) * n
        end
    end
    n_e_edge = max(n_pos_charge - n_neg_charge, eps(Float64))
    αs = n_neg_charge / n_e_edge

    # Electronegative correction factor from Ridenti et al (2025) Eq. (18)
    # Collapses to 1.0 when αs = 0 (no negative ions → classical Bohm)
    Te_eff_factor = (1 + αs) / (1 + γ * αs)

    # Neutral inlet density
    un = fluids.continuity[].const_velocity
    neutral_density = (mdot_a + ingestion_flow_rate) / cache.channel_area[1] / un

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
            0.0
        end
    else
        1.0
    end

    @inbounds for fluid in fluids.isothermal
        mi = fluid.species.element.m
        Z = fluid.species.Z

        interior_density = fluid.density[2]
        interior_flux = fluid.momentum[2]
        interior_velocity = interior_flux / interior_density

        if Z > 0
            # Electronegativity correction enters via Te_eff_factor.
            sound_speed = sqrt((kTi_J + Z * kTe_J * Te_eff_factor) / mi)
            boundary_velocity = -bohm_factor * sound_speed

            if interior_velocity <= -sound_speed
                # Supersonic outflow → pure Neumann
                boundary_density = interior_density
                boundary_flux = interior_flux
            else
                # Subsonic → drive to Bohm speed via Riemann invariants
                J⁻ = interior_velocity - sound_speed * log(interior_density)
                J⁺ = 2 * boundary_velocity - J⁻
                boundary_density = exp(0.5 * (J⁺ - J⁻) / sound_speed)
                boundary_flux = boundary_velocity * boundary_density
            end

            # send outflowing positive-ion flux back as neutrals
            neutral_density -= boundary_flux / un

        else
            # Negative ions: repelled by sheath, NOT beamed.
            sound_speed = sqrt(kTi_J / mi)
            boundary_velocity = 0.0

            if interior_velocity >= 0.0
                # Already flowing away from anode → Neumann
                boundary_density = interior_density
                boundary_flux = interior_density * interior_velocity
            else
                # Would leave domain → clamp to zero velocity
                J⁻ = interior_velocity - sound_speed * log(interior_density)
                J⁺ = -J⁻
                boundary_density = exp(0.5 * (J⁺ - J⁻) / sound_speed)
                boundary_flux = 0.0
            end
        end

        # Extrapolate to ghost cell
        fluid.density[1] = 2 * boundary_density - fluid.density[2]
        fluid.momentum[1] = 2 * boundary_flux - fluid.momentum[2]
    end

    # Extrapolate neutral density to ghost cell
    fluids.continuity[].density[1] = 2 * neutral_density - fluids.continuity[].density[2]

    return
end

function apply_right_boundary!(fluids)
    @inbounds for fluid in fluids.continuity
        fluid.density[end] = fluid.density[end - 1]
    end

    @inbounds for fluid in fluids.isothermal
        if fluid.species.Z > 0
            fluid.density[end] = fluid.density[end - 1]
            fluid.momentum[end] = fluid.momentum[end - 1]
        else
            interior_density = fluid.density[end - 1]
            interior_flux = fluid.momentum[end - 1]
            interior_velocity = interior_flux / interior_density

            if interior_velocity > 0
                fluid.density[end] = interior_density
                fluid.momentum[end] = interior_flux
            else
                fluid.density[end] = MIN_NUMBER_DENSITY * fluid.species.element.m
                fluid.momentum[end] = MIN_NUMBER_DENSITY * fluid.species.element.m * interior_velocity
            end
        end
    end

    return
end

#===============================================================================
Heavy species source terms
===============================================================================#

function apply_reactions!(fluid_arr, params)
    (;
        ei_reactions,
        ei_reactant_indices,
        ei_product_indices,
        cache, landmark,
    ) = params

    rxns = zip(
        ei_reactions, ei_reactant_indices, ei_product_indices,
    )

    return apply_reactions!(fluid_arr, rxns, cache, landmark)
end

function apply_reactions!(fluids, rxns, cache, landmark)
    (; inelastic_losses, νiz, ϵ, ne, K) = cache

    # Zero ionization frequency and inelastic losses and compute electron density
    @inbounds begin
        # Update electron density (TODO check if this is optimal)
        ne .= 0.0
        for fluid in fluids
            for i in eachindex(ne)
                ne[i] += fluid.species.Z * fluid.density[i] / fluid.species.element.m
            end
        end
        @. ne = max(ne, MIN_NUMBER_DENSITY)
        νiz .= 0.0
        inelastic_losses .= 0.0
        @. ϵ = cache.nϵ / cache.ne
        if !landmark
            @. ϵ += K
        end
    end

    dt_max = Inf
    for (rxn, reactant_index, product_index) in rxns
        # Temp storage for reaction calculations
        rxn_cache = (cache.cell_cache_1, cache.cell_cache_2)

        # Apply single reaction
        _dt = apply_reaction!(fluids, reactant_index, product_index, rxn.product_coeffs, rxn_cache, ne, ϵ, rxn, νiz, inelastic_losses, landmark)
        dt_max = min(_dt, dt_max)
    end

    cache.dt_iz[] = dt_max
    return
end

function apply_reaction!(fluids, reactant_index, product_index, product_coeffs, rxn_cache, ne, ϵ, rxn, νiz, inelastic_losses, landmark)
    dt_max = Inf
    reactant = fluids[reactant_index]
    reactant_velocity = reactant.const_velocity
    inv_m = 1 / reactant.species.element.m

    # Extract temp caches
    dens_cache, mom_cache = rxn_cache
    ncells = length(dens_cache)

    # Compute reaction rate and adjust reactant properties
    @inbounds @simd for i in 2:(ncells - 1)
        r = rate_coeff(rxn, ϵ[i])
        ρ_reactant = reactant.density[i]
        ρdot = reaction_rate(r, ne[i], ρ_reactant)
        dt_max = min(dt_max, ρ_reactant / ρdot)
        ndot = ρdot * inv_m
        νiz[i] += ndot / ne[i]
        inelastic_losses[i] += ndot * rxn.energy

        # Change in density due to ionization
        reactant.dens_ddt[i] -= ρdot

        # Store density changes in cache
        dens_cache[i] = ndot

        if !landmark
            if reactant.type != _ContinuityOnly
                # Momentum transfer due to ionization
                reactant_velocity = reactant.momentum[i] / ρ_reactant
                reactant.mom_ddt[i] -= ρdot * reactant_velocity
            end

            # Store momentum change in cache
            mom_cache[i] = ndot * reactant_velocity
        else
            mom_cache[i] = 0.0
        end
    end

    # Iterate products and add mass/momentum as needed
    @inbounds for (prod_ind, prod_coeff) in zip(product_index, product_coeffs)
        product = fluids[prod_ind]
        prod_mass = product.species.element.m

        @simd for i in 2:(ncells - 1)
            product.dens_ddt[i] += prod_mass * prod_coeff * dens_cache[i]
            product.mom_ddt[i] += prod_mass * prod_coeff * mom_cache[i]
        end
    end

    return dt_max
end

@inline reaction_rate(rate_coeff, ne, n_reactant) = rate_coeff * ne * n_reactant

function apply_ion_acceleration!(fluids::Vector{FluidContainer}, grid, cache)
    dt_max = Inf

    @inbounds for fluid in fluids
        Z = fluid.species.Z
        m = fluid.species.element.m
        qe_m = Z * e / m

        @simd for i in 2:(length(fluid.dens_ddt) - 1)
            qE_m = -qe_m * cache.∇ϕ[i]
            dz = grid.dz_cell[i]

            Q_accel = qE_m * fluid.density[i]
            dt_max = min(dt_max, abs(dz / qE_m))
            fluid.mom_ddt[i] += Q_accel
        end
    end

    return cache.dt_E[] = sqrt(dt_max)
end

function apply_ion_wall_losses!(fluid_containers, params)
    (; thruster, cache, grid, transition_length, wall_loss_scale) = params
    (; continuity, isothermal) = fluid_containers

    geometry = thruster.geometry
    L_ch = geometry.channel_length
    inv_Δr = inv(geometry.outer_radius - geometry.inner_radius)
    h = wall_loss_scale * edge_to_center_density_ratio()

    neutral_fluid = continuity[1]
    @inbounds for ion_fluid in isothermal
        m = ion_fluid.species.element.m
        Z = ion_fluid.species.Z
        qe_m = Z * e / m

        for i in 2:(length(ion_fluid.density) - 1)
            if Z > 0
                u_bohm = sqrt(qe_m * cache.Tev[i])
                in_channel = linear_transition(grid.cell_centers[i], L_ch, transition_length, 1.0, 0.0)
                νiw = in_channel * u_bohm * inv_Δr * h

                density_loss = ion_fluid.density[i] * νiw
                momentum_loss = ion_fluid.momentum[i] * νiw

                # Neutrals gain density due to ion recombination at the walls
                neutral_fluid.dens_ddt[i] += density_loss
                ion_fluid.dens_ddt[i] -= density_loss
                ion_fluid.mom_ddt[i] -= momentum_loss
            end
        end
    end

    return
end
