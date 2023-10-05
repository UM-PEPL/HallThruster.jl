function left_boundary_state!(bc_state, U, params)
    (;index, config,) = params
    (;Tev, channel_area) = params.cache
    mi = config.propellant.m
    Ti = config.ion_temperature

    un = config.neutral_velocity
    mdot_a = config.anode_mass_flow_rate

    if config.anode_boundary_condition == :sheath
        Vs = params.cache.Vs[]
        # Compute sheath potential
        electron_repelling_sheath = Vs > 0
        if electron_repelling_sheath
            # Ion attracting/electron-repelling sheath, ions in pre-sheath attain reduced Bohm speed
            Vs_norm = Vs / Tev[1]
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
    bc_state[index.ρn[1]] = mdot_a / channel_area[1] / un

    # Add ingested mass flow rate at anode
    if config.solve_background_neutrals
        bc_state[index.ρn[1]] += params.background_neutral_density * params.background_neutral_velocity / un
    end

    @inbounds for Z in 1:params.config.ncharge
        interior_density = U[index.ρi[Z],   begin+1]
        interior_flux    = U[index.ρiui[Z], begin+1]
        interior_velocity = interior_flux / interior_density

        sound_speed = sqrt((kB * Ti + Z * e * Tev[1]) / mi)  # Sound speed considering electron pressure-coupled terms
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

            #boundary_flux = interior_flux
            #boundary_density = boundary_flux / boundary_velocity
        end

        bc_state[index.ρn[1]] -= boundary_flux / un
        bc_state[index.ρi[Z]] = boundary_density
        bc_state[index.ρiui[Z]] = boundary_flux
    end
end

function right_boundary_state!(bc_state, U, params)
    (;index, fluids) = params

    @inbounds for j in 1:params.num_neutral_fluids
        # Use Neumann boundary conditions for all neutral fluids
        bc_state[index.ρn[j]] = U[index.ρn[j], end-1]
    end

    @inbounds for Z in 1:params.config.ncharge
        boundary_density = U[index.ρi[Z], end-1]
        bc_state[index.ρi[Z]]   = boundary_density        # Neumann BC for ion density at right boundary
        bc_state[index.ρiui[Z]] = U[index.ρiui[Z], end-1] # Neumann BC for ion flux at right boundary
    end
end
