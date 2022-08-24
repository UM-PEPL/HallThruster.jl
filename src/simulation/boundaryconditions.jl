function left_boundary_state!(bc_state, U, params)
    (;index, A_ch, config,) = params
    mi = config.propellant.m
    (;Tev) = params.cache

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
            χ = exp(-Vs_norm) / √(π * Vs_norm) / (1 + erf(sqrt(Vs_norm)))
            bohm_factor = inv(√(1 + χ))
        else
            # Ion-repelling sheath, ions have zero velocity at anode
            bohm_factor = 0.0
        end
    else
        bohm_factor = 1.0
    end

    # Precompute bohm velocity
    bohm_velocity = bohm_factor * sqrt(e * Tev[1] / mi)

    # Add inlet neutral density
    bc_state[index.ρn[1]] = mdot_a / A_ch / un

    if config.solve_background_neutrals
        # Background neutrals are accomodated and re-emitted as anode neutrals
        background_neutral_density = U[index.ρn[2], begin+1]
        background_neutral_velocity = params.background_neutral_velocity
        background_neutral_flux = -background_neutral_density * background_neutral_velocity
        bc_state[index.ρn[1]] += background_neutral_flux / un
        bc_state[index.ρn[2]] = U[index.ρn[2], begin+1]
    end

    @inbounds for Z in 1:params.config.ncharge
        boundary_density = U[index.ρi[Z],   begin+1]
        boundary_flux    = U[index.ρiui[Z], begin+1]
        boundary_velocity = boundary_flux / boundary_density

        # Enforce Bohm condition at left boundary
        boundary_velocity = min(-sqrt(Z) * bohm_velocity, boundary_velocity)

        recombination_density = -(boundary_density * boundary_velocity) / un

        bc_state[index.ρn[1]] += recombination_density
        bc_state[index.ρi[Z]] = boundary_density # Neumann BC for ion density at left boundary
        bc_state[index.ρiui[Z]] = boundary_velocity * boundary_density
    end
end

function right_boundary_state!(bc_state, U, params)
    (;index, fluids) = params

    @inbounds for j in 1:params.num_neutral_fluids
        # If neutral fluid is right-moving, use Neumann BC, otherwise Dirichlet BC
        if fluids[j].conservation_laws.u ≥ 0.0
            bc_state[index.ρn[j]] = U[index.ρn[j], end-1]
        else
            bc_state[index.ρn[j]] = params.background_neutral_density
        end
    end

    @inbounds for Z in 1:params.config.ncharge
        boundary_density = U[index.ρi[Z], end-1]
        bc_state[index.ρi[Z]]   = boundary_density        # Neumann BC for ion density at right boundary
        bc_state[index.ρiui[Z]] = U[index.ρiui[Z], end-1] # Neumann BC for ion flux at right boundary
    end
end
