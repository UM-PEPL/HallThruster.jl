function left_boundary_state!(bc_state, U, params)
    (;Te_L, index, A_ch, config, z_cell) = params
    mi = config.propellant.m
    (;Tev, ϕ) = params.cache

    c0, c1, c2 = second_deriv_coeffs(z_cell[1], z_cell[2], z_cell[3])

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

    if params.LANDMARK
        bohm_velocity = bohm_factor * sqrt(e * Te_L / mi)
    else
        bohm_velocity = bohm_factor * sqrt(e * Tev[1] / mi)
    end

    for Z in 1:params.config.ncharge
        boundary_density = U[index.ρi[Z], 2]
        boundary_flux = U[index.ρiui[Z], 2]
        boundary_velocity = boundary_flux / boundary_density

        ρ1, ρ2 = U[index.ρi[Z], 2], U[index.ρi[Z], 3]

        # Enforce Bohm condition
        boundary_velocity = min(-sqrt(Z) * bohm_velocity, boundary_velocity)

        #=if boundary_velocity < 0
            # Flux of density to the left boundary is conserved
            boundary_density = max(mi * params.config.min_number_density, boundary_flux / boundary_velocity)
        else
            # Second derivative of density is zero at the boundary
            boundary_density = -(c1 * ρ1 + c2 * ρ2) / c0
        end=#

        recombination_density = -(boundary_density * boundary_velocity) / un

        bc_state[index.ρn] += recombination_density
        bc_state[index.ρi[Z]] = boundary_density
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
