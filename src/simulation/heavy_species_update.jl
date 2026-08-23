function integrate_heavy_species!(fluid_containers, params, user_source::S, dt) where {S}
    # Strang-split the stiff, linear radiative subsystem around the SSPRK update.
    # Each half-step is an exact propagation of all coupled decay cascades.
    half_dt = 0.5 * dt
    apply_radiative_decay!(
        params.fluid_array,
        params.radiative_networks,
        half_dt,
    )
    # Do one timestep forward, returning `true` if we found a NaN or Inf
    step_heavy_species!(fluid_containers, params, user_source, dt) && return true
    apply_radiative_decay!(
        params.fluid_array,
        params.radiative_networks,
        half_dt,
    )
    limit_heavy_species!(fluid_containers)
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

As written, this requires three intermediate storage variables: k_{n1}, k_{n2}, and y_{n1}
We can reduce this to two using the following rearrangement

y_{n+1} = y_n / 2 + (y_n + h * k_{n1}) / 2 + h * k_{n2}
        = (y_n + y_{n1} + h * k_{n2}) / 2

With this, we do not need to store k_{n1} and k_{n2} separately and can instead reuse the same memory.

The implementation is split into two update kernels so each stage traverses every
fluid field only once. `update_first_stage!` saves y_n in the fluid caches while
forming the Euler predictor y_{n1}. After evaluating the derivatives at that
predictor, `update_second_stage!` combines the cached y_n, y_{n1}, and k_{n2} to
form the final state. Each kernel also performs its finite-value check during the
same traversal.
"""
function step_heavy_species!(fluid_containers, params, source::S, dt) where {S}
    # Evaluate k_{n1} at the initial state.
    compute_heavy_species_derivatives!(fluid_containers, params, source)
    params.cache.inelastic_losses_stage .= params.cache.inelastic_losses

    # Cache y_n, form the Euler predictor y_{n1}, and validate it in one pass.
    update_first_stage!(fluid_containers, dt) && return true
    # The second derivative evaluation must see a physically admissible predictor.
    limit_heavy_species!(fluid_containers)

    # Evaluate k_{n2} at y_{n1}, reusing the derivative arrays that held k_{n1}.
    compute_heavy_species_derivatives!(fluid_containers, params, source)

    # Match the energy sink to the same Heun-averaged reaction rate used for
    # the species source terms.
    @. params.cache.inelastic_losses = 0.5 * (
        params.cache.inelastic_losses_stage + params.cache.inelastic_losses
    )

    # Combine cached y_n with y_{n1} and k_{n2}, validating the result in one pass.
    update_second_stage!(fluid_containers, dt) && return true
    limit_heavy_species!(fluid_containers)

    return false
end

function update_first_stage!(fluid_containers, dt)
    invalid = false
    @inbounds for fluid in fluid_containers.continuity
        @simd for i in eachindex(fluid.density)
            density = fluid.density[i]
            new_density = density + dt * fluid.dens_ddt[i]
            fluid.dens_cache[i] = density
            fluid.density[i] = new_density
            invalid |= !isfinite(new_density)
        end
    end

    @inbounds for fluid in fluid_containers.isothermal
        @simd for i in eachindex(fluid.density)
            density = fluid.density[i]
            momentum = fluid.momentum[i]
            new_density = density + dt * fluid.dens_ddt[i]
            new_momentum = momentum + dt * fluid.mom_ddt[i]
            fluid.dens_cache[i] = density
            fluid.mom_cache[i] = momentum
            fluid.density[i] = new_density
            fluid.momentum[i] = new_momentum
            invalid |= !(isfinite(new_density) && isfinite(new_momentum))
        end
    end

    return invalid
end

function update_second_stage!(fluid_containers, dt)
    invalid = false

    @inbounds for fluid in fluid_containers.continuity
        @simd for i in eachindex(fluid.density)
            new_density = 0.5 * (
                fluid.density[i] + fluid.dens_cache[i] + dt * fluid.dens_ddt[i]
            )
            fluid.density[i] = new_density
            invalid |= !isfinite(new_density)
        end
    end

    @inbounds for fluid in fluid_containers.isothermal
        @simd for i in eachindex(fluid.density)
            new_density = 0.5 * (
                fluid.density[i] + fluid.dens_cache[i] + dt * fluid.dens_ddt[i]
            )
            new_momentum = 0.5 * (
                fluid.momentum[i] + fluid.mom_cache[i] + dt * fluid.mom_ddt[i]
            )
            fluid.density[i] = new_density
            fluid.momentum[i] = new_momentum
            invalid |= !(isfinite(new_density) && isfinite(new_momentum))
        end
    end

    return invalid
end

# Populate dens_ddt and mom_ddt for all fluid containers
function compute_heavy_species_derivatives!(
        fluid_containers, params, source_heavy_species::S,
    ) where {S}
    (; cache, grid, ion_wall_losses, reconstruct) = params

    update_convective_terms!(fluid_containers, grid, reconstruct, cache.dlnA_dz)
    source_heavy_species(fluid_containers, params)
    apply_reactions!(params.fluid_array, params)
    apply_mutual_neutralization!(params)
    apply_associative_detachment!(params)

    apply_ion_acceleration!(fluid_containers.isothermal, grid, cache)

    if ion_wall_losses
        apply_ion_wall_losses!(params)
    end

    # Update maximum allowable timestep
    CFL = params.simulation.CFL
    min_dt_u = Inf
    for fluid in params.fluid_array
        min_dt_u = min(min_dt_u, fluid.max_timestep[])
    end

    # The empirical 0.799 stability limit applies to chemistry. Transport and
    # acceleration can use a higher user-specified CFL independently.
    chemistry_CFL = min(CFL, 0.799)
    cache.dt[] = min(
        chemistry_CFL * cache.dt_iz[],
        sqrt(CFL) * cache.dt_E[],
        CFL * min_dt_u,
    )

    return
end

function limit_heavy_species!(fluid_containers)
    @inbounds for fluid in fluid_containers.continuity
        min_density = MIN_NUMBER_DENSITY * fluid.species.element.m
        @simd for i in eachindex(fluid.density)
            fluid.density[i] = max(fluid.density[i], min_density)
        end
    end

    @inbounds for fluid in fluid_containers.isothermal
        min_density = MIN_NUMBER_DENSITY * fluid.species.element.m
        @simd for i in eachindex(fluid.density)
            if fluid.density[i] < min_density
                fluid.density[i] = min_density
                fluid.momentum[i] = 0.0
            end
        end
    end
    return
end

function update_heavy_species!(params)
    (; cache, propellants, anode_bc, ingestion_flow_rates) = params

    # Apply left boundary conditions per-propellant
    for (i, (propellant, fluids)) in enumerate(zip(propellants, params.fluids_by_propellant))
        apply_left_boundary!(fluids, propellant, cache, anode_bc, ingestion_flow_rates[i], params.landmark)
    end

    # Apply right boundary conditions for all propellants
    apply_right_boundary!(params.fluid_containers)

    # Update ion variables as seen by electrons
    update_heavy_species_cache!(params.fluid_containers, cache, params.grid, params.landmark)

    return
end

function update_heavy_species_cache!(fluids, cache, grid, landmark)
    (; nn, ne, Z_eff, ji, ϵ, nϵ, K, m_eff, avg_ion_vel, avg_neutral_vel) = cache

    isempty(fluids.continuity) && throw(ArgumentError(
        "At least one neutral heavy species is required to update the plasma state."
    ))
    isempty(fluids.isothermal) && throw(ArgumentError(
        "At least one charged heavy species is required to update the plasma state."
    ))
    fallback_ion_mass = first(fluids.isothermal).species.element.m

    @inbounds @simd for i in eachindex(ne)
        ne[i] = 0.0
        ji[i] = 0.0
        m_eff[i] = 0.0
        Z_eff[i] = 0.0
        nn[i] = 0.0
        avg_ion_vel[i] = 0.0
        avg_neutral_vel[i] = 0.0
    end

    # Compute neutral number density, summed over all electronic states
    @inbounds for fluid in fluids.continuity
        inv_m = inv(fluid.species.element.m)
        neutral_velocity = fluid.const_velocity
        @simd for i in eachindex(fluid.density)
            number_density = fluid.density[i] * inv_m
            nn[i] += number_density
            avg_neutral_vel[i] += number_density * neutral_velocity
        end
    end


    # Update plasma quantities
    @inbounds for fluid in fluids.isothermal
        inv_m = inv(fluid.species.element.m)
        Z = fluid.species.Z

        @simd for i in eachindex(fluid.density)
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

    @inbounds @simd for i in eachindex(ne)
        neutral_density = nn[i]
        avg_neutral_vel[i] = neutral_density > 0 ?
            avg_neutral_vel[i] / neutral_density : 0.0
        ne[i] = max(ne[i], MIN_NUMBER_DENSITY)

        ion_density = Z_eff[i]
        if ion_density > 0
            inv_ion_density = inv(ion_density)
            avg_ion_vel[i] *= inv_ion_density
            m_eff[i] *= inv_ion_density
            Z_eff[i] = ne[i] * inv_ion_density
        else
            # A zero-density cell can occur in user initial conditions or a
            # restart before density limiting runs. Keep derived quantities
            # finite until the normal population floor is applied.
            avg_ion_vel[i] = 0.0
            m_eff[i] = fallback_ion_mass
            Z_eff[i] = 1.0
        end

        ϵ[i] = nϵ[i] / ne[i]
        if !landmark
            ϵ[i] += K[i]
        end
    end

    return
end

#===============================================================================
Boundary conditions
===============================================================================#

function apply_left_boundary!(fluids, propellant, cache, anode_bc, ingestion_flow_rate, landmark = false)
    Te_L = cache.Tev[2]                  # eV
    Ti = propellant.ion_temperature_K    # K
    mdot_a = propellant.flow_rate_kg_s

    kTe_J = e * Te_L      # electron energy in Joules
    kTi_J = kB * Ti       # ion thermal energy in Joules

    if !landmark
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
    else
        Te_eff_factor = 1.0
    end

    # Neutral inlet density. Anode flow feeds the ground state only.
    neutral_fluid = ground_neutral(fluids)
    un = neutral_fluid.const_velocity
    neutral_density = (mdot_a + ingestion_flow_rate) / cache.channel_area[1] / un

    Vs = 0.0
    bohm_factor = if anode_bc == :sheath
        Vs = cache.Vs[]
        # Compute sheath potential
        electron_repelling_sheath = Vs > 0
        if electron_repelling_sheath
            # Ion attracting/electron-repelling sheath, ions in pre-sheath attain reduced Bohm speed
            Vs_norm = (Vs / Te_L + 1.0e-6)
            # Compute correction factor (see Hara, PSST 28 (2019))
            χ = exp(-Vs_norm) / √(π * Vs_norm) / (1 + myerf(sqrt(Vs_norm)))
            inv(√(1 + χ))
        else
            1.0
        end
    else
        1.0
    end

    @inbounds for fluid in fluids.isothermal
        mi = fluid.species.element.m
        Z = fluid.species.Z

        interior_density = fluid.density[2]
        interior_flux = fluid.momentum[2]
        interior_velocity = primitive_velocity(interior_flux, interior_density)

        if Z > 0
            # Electronegativity correction enters via Te_eff_factor.
            sound_speed = sqrt((kTi_J + Z * kTe_J * Te_eff_factor) / mi)
            boundary_velocity = -bohm_factor * sound_speed
            interior_density_safe = max(interior_density, MIN_NUMBER_DENSITY * mi)

            if interior_velocity <= -sound_speed
                # Supersonic outflow → pure Neumann
                boundary_density = interior_density
                boundary_flux = interior_flux
            else
                # Subsonic outflow, need to drive the flow toward sonic
                # For the isothermal Euler equations, the Riemann invariants are
                # J⁺ = u + c ln ρ
                # J⁻ = u - c ln ρ
                # For the boundary condition, we take c = u_bohm and use J⁻ to set the boundary density.

                # J⁻ from interior (outgoing)
                J⁻ = interior_velocity - sound_speed * log(interior_density_safe)

                # Set boundary velocity to Bohm, use J⁻ to get boundary density
                # J⁻ = boundary_velocity - sound_speed * log(boundary_density)
                # → log(boundary_density) = (boundary_velocity - J⁻) / sound_speed
                boundary_density = exp((boundary_velocity - J⁻) / sound_speed)
                boundary_flux = boundary_velocity * boundary_density
            end

            # send outflowing positive-ion flux back as ground-state neutrals
            neutral_density -= boundary_flux / un

        else
            KE = 0.5 * mi * interior_velocity^2
            barrier = abs(Z) * e * Vs

            if KE >= barrier
                # Ion reaches anode with reduced speed
                KE_boundary = KE - barrier
                boundary_velocity = -sqrt(2 * KE_boundary / mi)  # toward anode
                boundary_density = interior_density
                boundary_flux = boundary_density * boundary_velocity
            else
                # Reflected
                boundary_density = interior_density
                boundary_flux = 0.0
            end
        end

        fluid.density[1] = boundary_density
        fluid.momentum[1] = boundary_flux
    end

    nm = neutral_fluid.species.element.m
    neutral_fluid.density[1] = max(neutral_density, MIN_NUMBER_DENSITY * nm)

    # Excited states have no anode inflow, so clamp their ghost cells to the floor
    @inbounds for fluid in fluids.continuity
        is_excited(fluid.species) || continue
        fluid.density[1] = MIN_NUMBER_DENSITY * fluid.species.element.m
    end

    return
end

function apply_right_boundary!(fluids)
    @inbounds for fluid in fluids.continuity
        fluid.density[end] = fluid.density[end - 1]
    end

    @inbounds for fluid in fluids.isothermal
        interior_density = fluid.density[end - 1]
        interior_flux = fluid.momentum[end - 1]
        interior_velocity = primitive_velocity(interior_flux, interior_density)
        mi = fluid.species.element.m

        if interior_velocity >= 0
            # Normal supersonic outflow — Neumann
            fluid.density[end] = interior_density
            fluid.momentum[end] = interior_flux
        else
            # inflow clamp
            fluid.density[end] = MIN_NUMBER_DENSITY * mi
            fluid.momentum[end] = MIN_NUMBER_DENSITY * mi * interior_velocity
        end
    end

    return
end

#===============================================================================
Heavy species source terms
===============================================================================#

"""Precomputed metadata for one electron-impact reaction channel."""
struct ElectronImpactChannel
    reaction::ElectronImpactReaction
    product_indices::Vector{Int}
    product_mass_ratios::Vector{Float64}
    is_ionizing::Bool
    is_excitation::Bool
end

"""
Electron-impact channels sharing one reactant and loss-frequency accumulator.
"""
struct ElectronImpactGroup
    reactant_index::Int
    inverse_reactant_mass::Float64
    carries_momentum::Bool
    channels::Vector{ElectronImpactChannel}
end

function apply_reactions!(fluid_arr, params)
    return apply_reaction_groups!(
        fluid_arr, params.reaction_groups, params.cache, params.landmark,
    )
end

"""Return the shared lookup-table index limit, or -1 when tables differ."""
function common_rate_index_limit(groups)
    limit = -1
    for group in groups
        for channel in group.channels
            channel_limit = length(channel.reaction.rate_coeffs) - 2
            channel_limit >= 0 || return -1
            if limit < 0
                limit = channel_limit
            elseif channel_limit != limit
                return -1
            end
        end
    end
    return limit
end

function prepare_reaction_state!(fluids, cache, landmark, groups)
    (; inelastic_losses, νiz, νex_explicit, ϵ, ne, K) = cache

    # Recompute electron density for the current RK stage. Neutral fluids do not
    # contribute, which becomes increasingly useful with multiple propellants.
    fill!(ne, 0.0)
    @inbounds for fluid in fluids
        Z = fluid.species.Z
        iszero(Z) && continue
        charge_to_mass = Z / fluid.species.element.m
        @simd for i in eachindex(ne)
            ne[i] += charge_to_mass * fluid.density[i]
        end
    end

    # Initialize reaction outputs and electron energy in the same traversal.
    @inbounds @simd for i in eachindex(ne)
        electron_density = max(ne[i], MIN_NUMBER_DENSITY)
        ne[i] = electron_density
        νiz[i] = 0.0
        νex_explicit[i] = 0.0
        inelastic_losses[i] = 0.0
        ϵ[i] = cache.nϵ[i] / electron_density + (landmark ? 0.0 : K[i])
    end
    # When lookup tables share their unit-spaced coordinate, clamp the index and
    # compute its interpolation fraction once per cell instead of once per table.
    rate_index_limit = if hasproperty(cache, :reaction_rate_index_limit)
        cache.reaction_rate_index_limit[]
    else
        common_rate_index_limit(groups)
    end
    has_lookup_cache = hasproperty(cache, :reaction_rate_indices) &&
        hasproperty(cache, :reaction_rate_fractions)
    if rate_index_limit >= 0 && has_lookup_cache
        indices = cache.reaction_rate_indices
        fractions = cache.reaction_rate_fractions
        @inbounds @simd for i in eachindex(ϵ)
            energy = ϵ[i]
            if isfinite(energy)
                index = clamp(Base.unsafe_trunc(Int, energy), 0, rate_index_limit)
                indices[i] = index
                fractions[i] = energy - index
            else
                indices[i] = 0
                fractions[i] = 0.0
            end
        end
        return true
    end
    return false
end

# Electronic excitation preserves the gas and charge state while changing its
# explicitly tracked level. Other charge-conserving reactions may dissociate.
@inline function _is_electronic_excitation(rxn)
    length(rxn.products) == 1 || return false
    only(rxn.product_coeffs) == 1 || return false
    product = only(rxn.products)
    return product.element.formula == rxn.reactant.element.formula &&
        product.Z == rxn.reactant.Z &&
        product.excited_level != rxn.reactant.excited_level
end

# A reaction is ionizing if any product's charge state differs from the reactant's.
# Excitation and charge-conserving dissociation are not, and must not enter νiz.
@inline function _is_ionizing(fluids, reactant_index, product_index)
    reactant_Z = fluids[reactant_index].species.Z
    for prod_ind in product_index
        if fluids[prod_ind].species.Z != reactant_Z
            return true
        end
    end
    return false
end

function build_electron_impact_groups(
        reactions, reactant_indices, product_indices, fluids,
    )
    grouped_channels = OrderedDict{Int, Vector{ElectronImpactChannel}}()
    for (reaction, reactant_index, products) in zip(
            reactions, reactant_indices, product_indices,
        )
        reactant_mass = fluids[reactant_index].species.element.m
        product_mass_ratios = [
            fluids[product_index].species.element.m * coefficient / reactant_mass
                for (product_index, coefficient) in zip(products, reaction.product_coeffs)
        ]
        channel = ElectronImpactChannel(
            reaction,
            products,
            product_mass_ratios,
            _is_ionizing(fluids, reactant_index, products),
            _is_electronic_excitation(reaction),
        )
        push!(get!(grouped_channels, reactant_index, ElectronImpactChannel[]), channel)
    end

    return [
        ElectronImpactGroup(
            reactant_index,
            inv(fluids[reactant_index].species.element.m),
            fluids[reactant_index].type != _ContinuityOnly,
            channels,
        ) for (reactant_index, channels) in pairs(grouped_channels)
    ]
end

function apply_reaction_groups!(fluids, groups, cache, landmark)
    use_cached_coordinates = prepare_reaction_state!(
        fluids, cache, landmark, groups,
    )
    reaction_rate_indices = use_cached_coordinates ? cache.reaction_rate_indices : nothing
    reaction_rate_fractions = use_cached_coordinates ? cache.reaction_rate_fractions : nothing
    loss_frequency = cache.reaction_loss_frequency
    max_loss_frequency = 0.0

    for group in groups
        isempty(group.channels) && continue
        last_channel = length(group.channels)
        if last_channel == 1
            group_max = apply_reaction_channel!(
                fluids, group, first(group.channels), cache, landmark,
                loss_frequency, reaction_rate_indices, reaction_rate_fractions,
                Val(true), Val(true),
            )
        else
            apply_reaction_channel!(
                fluids, group, first(group.channels), cache, landmark,
                loss_frequency, reaction_rate_indices, reaction_rate_fractions,
                Val(true), Val(false),
            )
            for channel_index in 2:(last_channel - 1)
                apply_reaction_channel!(
                    fluids, group, group.channels[channel_index], cache, landmark,
                    loss_frequency, reaction_rate_indices, reaction_rate_fractions,
                    Val(false), Val(false),
                )
            end
            group_max = apply_reaction_channel!(
                fluids, group, last(group.channels), cache, landmark,
                loss_frequency, reaction_rate_indices, reaction_rate_fractions,
                Val(false), Val(true),
            )
        end
        max_loss_frequency = max(max_loss_frequency, group_max)
    end

    cache.dt_iz[] = max_loss_frequency > 0 ? inv(max_loss_frequency) : Inf
    return nothing
end

function apply_reaction_channel!(
        fluids, group, channel, cache, landmark,
        loss_frequency, reaction_rate_indices, reaction_rate_fractions,
        ::Val{FIRST}, ::Val{LAST},
    ) where {FIRST, LAST}
    (; inelastic_losses, νiz, νex_explicit, ϵ, ne) = cache
    reaction = channel.reaction
    reactant = fluids[group.reactant_index]
    density_loss_cache = cache.cell_cache_1
    ncells = length(density_loss_cache)
    group_max = 0.0

    @inbounds @simd for cell in 2:(ncells - 1)
        rate = if isnothing(reaction_rate_indices)
            rate_coeff(reaction, ϵ[cell])
        else
            cached_rate_coeff(
                reaction, reaction_rate_indices[cell], reaction_rate_fractions[cell],
            )
        end
        reactant_density = reactant.density[cell]
        destruction_frequency = rate * ne[cell]
        density_loss = destruction_frequency * reactant_density
        reaction_frequency = rate * reactant_density * group.inverse_reactant_mass
        reactant_velocity = if landmark
            0.0
        elseif group.carries_momentum
            reactant.vel_prim[cell]
        else
            reactant.const_velocity
        end

        positive_loss_frequency = density_loss > 0 ? destruction_frequency : 0.0
        if FIRST
            loss_frequency[cell] = positive_loss_frequency
        else
            loss_frequency[cell] += positive_loss_frequency
        end
        channel.is_ionizing && (νiz[cell] += reaction_frequency)
        channel.is_excitation && (νex_explicit[cell] += reaction_frequency)
        inelastic_losses[cell] +=
            density_loss * group.inverse_reactant_mass * reaction.energy
        reactant.dens_ddt[cell] -= density_loss

        if !landmark
            if group.carries_momentum
                reactant.mom_ddt[cell] -= density_loss * reactant_velocity
            end
        end
        density_loss_cache[cell] = density_loss

        LAST && (group_max = max(group_max, loss_frequency[cell]))
    end

    @inbounds for (product_index, mass_ratio) in zip(
            channel.product_indices, channel.product_mass_ratios,
        )
        product = fluids[product_index]
        if landmark
            @simd for cell in 2:(ncells - 1)
                product.dens_ddt[cell] += mass_ratio * density_loss_cache[cell]
            end
        elseif group.carries_momentum
            reactant_velocity_cache = reactant.vel_prim
            @simd for cell in 2:(ncells - 1)
                mass_source = mass_ratio * density_loss_cache[cell]
                product.dens_ddt[cell] += mass_source
                product.mom_ddt[cell] += mass_source * reactant_velocity_cache[cell]
            end
        else
            @simd for cell in 2:(ncells - 1)
                mass_source = mass_ratio * density_loss_cache[cell]
                product.dens_ddt[cell] += mass_source
                product.mom_ddt[cell] += mass_source * reactant.const_velocity
            end
        end
    end
    return group_max
end

function apply_reaction!(
        fluids, reactant_index, product_index, product_coeffs, rxn_cache,
        ne, ϵ, rxn, νiz, νex_explicit, inelastic_losses, landmark,
        loss_frequency = nothing, reaction_rate_indices = nothing,
    )
    max_destruction_frequency = 0.0
    reactant = fluids[reactant_index]
    reactant_velocity = reactant.const_velocity
    inv_m = 1 / reactant.species.element.m

    # Only ionizing channels contribute to νiz; all channels contribute inelastic losses
    is_ionizing = _is_ionizing(fluids, reactant_index, product_index)
    is_excitation = _is_electronic_excitation(rxn)

    # Extract temp caches
    dens_cache, mom_cache = rxn_cache
    ncells = length(dens_cache)

    # Compute reaction rate and adjust reactant properties
    @inbounds @simd for i in 2:(ncells - 1)
        r = if isnothing(reaction_rate_indices)
            rate_coeff(rxn, ϵ[i])
        else
            rate_coeff(rxn, ϵ[i], reaction_rate_indices[i])
        end
        ρ_reactant = reactant.density[i]
        destruction_frequency = r * ne[i]
        ρdot = destruction_frequency * ρ_reactant
        ndot = ρdot * inv_m
        if ρdot > 0
            if isnothing(loss_frequency)
                max_destruction_frequency = max(
                    max_destruction_frequency, destruction_frequency,
                )
            else
                loss_frequency[i] += destruction_frequency
            end
        end
        reaction_frequency = r * ρ_reactant * inv_m
        if is_ionizing
            νiz[i] += reaction_frequency
        end
        if is_excitation
            νex_explicit[i] += reaction_frequency
        end
        inelastic_losses[i] += ndot * rxn.energy

        # Change in density due to this reaction
        reactant.dens_ddt[i] -= ρdot

        # Store density changes in cache
        dens_cache[i] = ndot

        if !landmark
            if reactant.type != _ContinuityOnly
                # Momentum transfer due to ionization
                reactant_velocity = primitive_velocity(reactant.momentum[i], ρ_reactant)
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

    return inv(max_destruction_frequency)
end

@inline reaction_rate(rate_coeff, ne, n_reactant) = rate_coeff * ne * n_reactant

#===============================================================================
Mutual neutralization
===============================================================================#

# Ion-ion neutralization rate coefficient (lit. ~4e-14 to ~1e-12 m^3/s).
const _K_MUTUAL_NEUTRALIZATION = 1.0e-12   # m^3/s

# Aggregated neutral-induced detachment rate (A^- + N -> A + N + e^-); n_e-independent, dominates the cold shoulder where MN is weak.
const _K_NEUTRAL_DETACHMENT = 1.0e-12   # m^3/s

# Apply mutual neutralization A^- + B^+ -> neutral A + neutral B: the dominant n_e-independent negative-ion sink that bounds α = n_-/n_e; added as a direct source term since the electron-impact framework can't express it.
function apply_mutual_neutralization!(params)
    k_MN = _K_MUTUAL_NEUTRALIZATION
    propellant_groups = params.fluids_by_propellant
    cache = params.cache

    dt_max = Inf

    @inbounds for neg_group in propellant_groups
        # Ground-state neutral that receives the neutralized negative ion (e.g. H^- -> H).
        neutral_neg = ground_neutral(neg_group)
        for neg_fluid in neg_group.isothermal
            neg_fluid.species.Z < 0 || continue
            m_neg = neg_fluid.species.element.m
            inv_m_neg = inv(m_neg)

            for pos_group in propellant_groups
                # Ground-state neutral that receives the neutralized positive ion (e.g. H2O^+ -> H2O).
                neutral_pos = ground_neutral(pos_group)
                for pos_fluid in pos_group.isothermal
                    pos_fluid.species.Z > 0 || continue
                    m_pos = pos_fluid.species.element.m
                    inv_m_pos = inv(m_pos)

                    for i in eachindex(neg_fluid.dens_ddt)
                        ρ_neg = neg_fluid.density[i]
                        ρ_pos = pos_fluid.density[i]
                        n_neg = ρ_neg * inv_m_neg
                        n_pos = ρ_pos * inv_m_pos

                        R = k_MN * n_neg * n_pos      # 1/m^3/s

                        ρdot_neg = R * m_neg
                        ρdot_pos = R * m_pos

                        # Ions destroyed; parent neutrals gain the mass.
                        neg_fluid.dens_ddt[i] -= ρdot_neg
                        pos_fluid.dens_ddt[i] -= ρdot_pos
                        neutral_neg.dens_ddt[i] += ρdot_neg
                        neutral_pos.dens_ddt[i] += ρdot_pos

                        # Destroyed ions carry off their drift momentum; products are near-thermal.
                        u_neg = neg_fluid.momentum[i] / max(ρ_neg, eps())
                        u_pos = pos_fluid.momentum[i] / max(ρ_pos, eps())
                        neg_fluid.mom_ddt[i] -= ρdot_neg * u_neg
                        pos_fluid.mom_ddt[i] -= ρdot_pos * u_pos

                        # CFL: dt < density / |destruction rate| per species.
                        if ρdot_neg > eps()
                            dt_max = min(dt_max, ρ_neg / ρdot_neg)
                        end
                        if ρdot_pos > eps()
                            dt_max = min(dt_max, ρ_pos / ρdot_pos)
                        end
                    end
                end
            end
        end
    end

    # Tighten the chemistry CFL with the MN destruction-rate constraint.
    cache.dt_iz[] = min(cache.dt_iz[], dt_max)

    return nothing
end

#===============================================================================
Associative / collisional detachment by neutrals
===============================================================================#

#Apply neutral-induced detachment A^- + N -> A + N + e^- at rate k_AD·n(A^-)·n_neutral (from cache.nn): n_e-independent, so it dominates the cold shoulder where MN is weak; the freed electron is implicit in n_e = n_+ - n_-.

function apply_associative_detachment!(params)
    k_AD = _K_NEUTRAL_DETACHMENT
    propellant_groups = params.fluids_by_propellant
    cache = params.cache

    dt_max = Inf

    @inbounds for neg_group in propellant_groups
        neutral_neg = ground_neutral(neg_group)
        for neg_fluid in neg_group.isothermal
            neg_fluid.species.Z < 0 || continue
            m_neg = neg_fluid.species.element.m
            inv_m_neg = inv(m_neg)

            for i in eachindex(neg_fluid.dens_ddt)
                ρ_neg = neg_fluid.density[i]
                n_neg = ρ_neg * inv_m_neg
                n_n = cache.nn[i]            # total neutral number density

                R = k_AD * n_neg * n_n       # 1/m^3/s
                ρdot = R * m_neg

                # H^- destroyed; the H atom joins the H neutral fluid (partner neutral unchanged here).
                neg_fluid.dens_ddt[i] -= ρdot
                neutral_neg.dens_ddt[i] += ρdot

                u_neg = neg_fluid.momentum[i] / max(ρ_neg, eps())
                neg_fluid.mom_ddt[i] -= ρdot * u_neg

                if ρdot > eps()
                    dt_max = min(dt_max, ρ_neg / ρdot)
                end
            end
        end
    end

    cache.dt_iz[] = min(cache.dt_iz[], dt_max)

    return nothing
end

function apply_ion_acceleration!(fluids::Vector{FluidContainer}, grid, cache)
    max_abs_qe_m = 0.0

    @inbounds for fluid in fluids
        Z = fluid.species.Z
        m = fluid.species.element.m
        qe_m = Z * e / m
        max_abs_qe_m = max(max_abs_qe_m, abs(qe_m))

        @simd for i in 2:(length(fluid.dens_ddt) - 1)
            qE_m = -qe_m * cache.∇ϕ[i]
            fluid.mom_ddt[i] += qE_m * fluid.density[i]
        end
    end

    # The most restrictive acceleration timestep comes from the ion with the
    # largest |q/m|, so reduce over the grid once rather than once per species.
    dt_max = Inf
    @inbounds @simd for i in 2:(length(grid.dz_cell) - 1)
        qE_m = max_abs_qe_m * cache.∇ϕ[i]
        # Skip non-finite ∇ϕ so NaN cannot poison dt_E through the reduction.
        if isfinite(qE_m)
            dt_max = min(dt_max, abs(grid.dz_cell[i] / qE_m))
        end
    end

    return cache.dt_E[] = isfinite(dt_max) ? sqrt(dt_max) : Inf
end

function prepare_ion_wall_losses!(params)
    (; thruster, cache, wall_loss_scale) = params
    geometry = thruster.geometry
    inv_Δr = inv(geometry.outer_radius - geometry.inner_radius)
    h = wall_loss_scale * edge_to_center_density_ratio()
    wall_loss_base = cache.cell_cache_1
    wall_cells = 2:params.last_wall_cell

    # This cell-dependent part is common to every ion species: it contains the
    # local electron temperature, wall transition, geometry, and sheath-density
    # correction. Compute it once per derivative evaluation and reuse the cache
    # while applying the species-dependent charge-to-mass scaling below.
    @inbounds @simd for i in wall_cells
        wall_loss_base[i] =
            cache.wall_transition[i] * sqrt(e * cache.Tev[i]) * inv_Δr * h
    end

    return wall_loss_base, wall_cells
end

function apply_ion_wall_losses!(params)
    # Preparing outside the propellant loop avoids repeating the shared cell work
    # for molecular propellants with many ion fluids or reaction products.
    wall_loss_base, wall_cells = prepare_ion_wall_losses!(params)
    for fluids in params.fluids_by_propellant
        apply_ion_wall_losses!(fluids, wall_loss_base, wall_cells)
    end
    return
end

function apply_ion_wall_losses!(fluid_containers, params)
    # Retain the single-container entry point while using the same split between
    # shared cell work and species-specific losses.
    wall_loss_base, wall_cells = prepare_ion_wall_losses!(params)
    return apply_ion_wall_losses!(fluid_containers, wall_loss_base, wall_cells)
end

function apply_ion_wall_losses!(fluid_containers, wall_loss_base, wall_cells)
    (; continuity, isothermal) = fluid_containers

    neutral_fluid = ground_neutral(fluid_containers)
    @inbounds for ion_fluid in isothermal
        Z = ion_fluid.species.Z

        # Do not apply wall losses to negative ions, as they are repelled from the positive pre-sheath.
        if Z < 0
            continue
        end

        # Complete the ion wall-loss frequency by applying sqrt(Z / m) to the
        # shared cell-dependent factor prepared above.
        species_scale = sqrt(Z / ion_fluid.species.element.m)

        for i in wall_cells
            νiw = wall_loss_base[i] * species_scale

            density_loss = ion_fluid.density[i] * νiw
            momentum_loss = ion_fluid.momentum[i] * νiw

            # Neutrals gain density due to ion recombination at the walls
            neutral_fluid.dens_ddt[i] += density_loss
            ion_fluid.dens_ddt[i] -= density_loss
            ion_fluid.mom_ddt[i] -= momentum_loss
        end
    end

    return
end
