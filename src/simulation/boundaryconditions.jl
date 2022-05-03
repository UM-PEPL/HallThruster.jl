function left_boundary_state!(bc_state, U, params)
    (;Te_L, index, A_ch, config, z_cell) = params
    mi = config.propellant.m
    (;Tev, ϕ) = params.cache

    un = config.neutral_velocity
    mdot_a = config.anode_mass_flow_rate
    bc_state[index.ρn] = mdot_a / A_ch / un

    ne = 0.0
    if config.LANDMARK
        bohm_factor = 1.0
    else
        Vs = params.ϕ_L - ϕ[1]
        if Vs < 0
            # Ion attracting/electron-repelling sheath, ions in pre-sheath attain reduced Bohm speed
            Vs_norm = Vs / Tev[1]
            # Compute correction factor (see Hara, PSST 28 (2019))
            χ = exp(Vs_norm) / √(-π * Vs_norm) / (1 + erf(sqrt(-Vs_norm)))
            bohm_factor = inv(√(1 + χ))
        else
            # Ion-repelling sheath, ions have zero velocity at anode
            bohm_factor = 0.0
        end
    end

    if params.config.LANDMARK
        bohm_velocity = bohm_factor * sqrt(e * Te_L / mi)
    else
        bohm_velocity = bohm_factor * sqrt(e * Tev[1] / mi)
    end

    for Z in 1:params.config.ncharge
        boundary_density = U[index.ρi[Z], 2]
        boundary_flux = U[index.ρiui[Z], 2]
        boundary_velocity = boundary_flux / boundary_density

        # Enforce Bohm condition
        boundary_velocity = min(-sqrt(Z) * bohm_velocity, boundary_velocity)

        recombination_density = -(boundary_density * boundary_velocity) / un

        bc_state[index.ρn] += recombination_density
        bc_state[index.ρi[Z]] = boundary_density # Neumann BC for ion density at left boundary
        bc_state[index.ρiui[Z]] = boundary_velocity * boundary_density

        ne += Z * boundary_density / mi
    end
end

function right_boundary_state!(bc_state, U, params)
    (;Te_R, index) = params
    mi = params.config.propellant.m

    bc_state[index.ρn] = U[index.ρn, end-1]

    ne = 0.0
    for Z in 1:params.config.ncharge
        boundary_density = U[index.ρi[Z], end-1]
        boundary_flux = U[index.ρiui[Z], end-1]
        bc_state[index.ρi[Z]] = boundary_density
        bc_state[index.ρiui[Z]] = boundary_flux

        ne += Z * boundary_density / mi
    end
end
