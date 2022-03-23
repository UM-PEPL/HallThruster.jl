function left_boundary_state!(bc_state, U, params)
    (;Te_L, index, mdot_a, A_ch, fluids) = params
    mi = params.config.propellant.m

    un = fluids[1].conservation_laws.u
    bc_state[index.ρn] = mdot_a / A_ch / un

    ne = 0.0
    bohm_velocity = sqrt(e * 2/3 * Te_L / mi)
    for Z in 1:params.config.ncharge
        boundary_density = U[index.ρi[Z]]
        boundary_flux = U[index.ρiui[Z]]
        boundary_velocity = boundary_flux / boundary_density

        # Enforce Bohm condition
        boundary_velocity = min(-sqrt(Z) * bohm_velocity, boundary_velocity)
        boundary_density = max(mi * params.config.min_number_density, boundary_flux / boundary_velocity)

        bc_state[index.ρn] -= boundary_flux / un
        bc_state[index.ρi[Z]] = boundary_density
        bc_state[index.ρiui[Z]] = boundary_flux

        ne += Z * boundary_density / mi
    end

    bc_state[index.nϵ] = ne * Te_L
end

function right_boundary_state!(bc_state, U, params)
    (;Te_R, index) = params
    mi = params.config.propellant.m

    bc_state[index.ρn] = U[index.ρn]

    ne = 0.0
    for Z in 1:params.config.ncharge
        boundary_density = U[index.ρi[Z], 1]
        boundary_flux = U[index.ρiui[Z]]
        bc_state[index.ρi[Z]] = boundary_density
        bc_state[index.ρiui[Z]] = boundary_flux

        ne += Z * boundary_density / mi
    end

    bc_state[index.nϵ] = ne * Te_R
end
