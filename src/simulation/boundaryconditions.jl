function left_boundary_state!(bc_state, U, params)
    index = params.index
    ncharge = params.ncharge
    mi = params.mi
    Ti = params.ion_temperature_K
    un = params.neutral_velocity
    mdot_a = params.anode_mass_flow_rate
    nn_B = params.background_neutral_density
    un_B = params.background_neutral_velocity
    neutral_ingestion_multiplier = params.neutral_ingestion_multiplier
    anode_bc = params.anode_bc

    ingestion_density = nn_B * un_B / un * neutral_ingestion_multiplier

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

    return @inbounds for Z in 1:ncharge
        boundary_density = U[index.ρi[Z], end - 1]
        bc_state[index.ρi[Z]] = boundary_density        # Neumann BC for ion density at right boundary
        bc_state[index.ρiui[Z]] = U[index.ρiui[Z], end - 1] # Neumann BC for ion flux at right boundary
    end
end
