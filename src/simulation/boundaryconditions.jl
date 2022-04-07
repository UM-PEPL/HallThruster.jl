function left_boundary_state!(bc_state, U, params)
    (;Te_L, index, A_ch, config, z_cell) = params
    mi = config.propellant.m

    c0, c1, c2 = second_deriv_coeffs(z_cell[1], z_cell[2], z_cell[3])

    un = config.neutral_velocity
    mdot_a = config.anode_mass_flow_rate
    bc_state[index.ρn] = mdot_a / A_ch / un

    ne = 0.0
    bohm_velocity = sqrt(e * Te_L / mi)
    for Z in 1:params.config.ncharge
        boundary_density = U[index.ρi[Z], 2]
        boundary_flux = U[index.ρiui[Z], 2]
        boundary_velocity = boundary_flux / boundary_density

        ρ1, ρ2 = U[index.ρi[Z], 2], U[index.ρi[Z], 3]

        # Enforce Bohm condition
        boundary_velocity = min(-sqrt(Z) * bohm_velocity, boundary_velocity)
        boundary_density = -(c1 * ρ1 + c2 * ρ2) / c0 #max(mi * params.config.min_number_density, boundary_flux / boundary_velocity)

        bc_state[index.ρn] -= (boundary_density * boundary_velocity) / un
        bc_state[index.ρi[Z]] = boundary_density
        bc_state[index.ρiui[Z]] = boundary_velocity * boundary_density

        ne += Z * boundary_density / mi
    end

    bc_state[index.nϵ] = ne * Te_L
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

    bc_state[index.nϵ] = ne * Te_R
end
