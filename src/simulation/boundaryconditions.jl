function left_boundary_state!(bc_state, U, params)
    (;index, A_ch, config) = params
    mi = config.propellant.m
    (;Tev) = params.cache

    un = config.neutral_velocity
    mdot_a = config.anode_mass_flow_rate

    if config.LANDMARK
        bohm_factor = 1.0
    else
        Vs = params.cache.Vs[]
        # Compute sheath potential
        electron_repelling_sheath = Vs > 0
        if electron_repelling_sheath
            # Ion attracting/electron-repelling sheath, ions in pre-sheath attain reduced Bohm speed
            Vs_norm = Vs / Tev[1]
            # Compute correction factor (see Hara, PSST 28 (2019))
            χ = exp(-Vs_norm) / √(π * Vs_norm) / (1 + erf(sqrt(Vs_norm)))
            bohm_factor = inv(√(1 + χ))
        else
            # Ion-repelling sheath, ions have zero velocity at anode
            bohm_factor = 0.0
        end
    end

    # Precompute bohm velocity
    bohm_velocity = bohm_factor * sqrt(e * Tev[1] / mi)

    # Add inlet neutral density
    bc_state[index.ρn] = mdot_a / A_ch / un

    @inbounds for Z in 1:params.config.ncharge
        boundary_density = U[index.ρi[Z],   begin+1]
        boundary_flux    = U[index.ρiui[Z], begin+1]
        boundary_velocity = boundary_flux / boundary_density

        # Enforce Bohm condition at left boundary
        boundary_velocity = min(-sqrt(Z) * bohm_velocity, boundary_velocity)

        recombination_density = -(boundary_density * boundary_velocity) / un

        bc_state[index.ρn] += recombination_density
        bc_state[index.ρi[Z]] = boundary_density # Neumann BC for ion density at left boundary
        bc_state[index.ρiui[Z]] = boundary_velocity * boundary_density
    end
end

function right_boundary_state!(bc_state, U, params)
    (;index) = params
    bc_state[index.ρn] = U[index.ρn, end-1]

    @inbounds for Z in 1:params.config.ncharge
        boundary_density = U[index.ρi[Z], end-1]
        bc_state[index.ρi[Z]]   = boundary_density        # Neumann BC for ion density at right boundary
        bc_state[index.ρiui[Z]] = U[index.ρiui[Z], end-1] # Neumann BC for ion flux at right boundary
    end
end
