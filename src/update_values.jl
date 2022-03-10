
function update_values!(integrator)
    (nvars, ncells) = size(integrator.u)
    @inbounds for j in 1:ncells, i in 1:nvars
        if isnan(integrator.u[i, j])
            println("NaN detected in variable $i in cell $j at time $(integrator.t)")
            terminate!(integrator, :NaNDetected)
        elseif isinf(integrator.u[i, j])
            println("Inf detected in variable $i in cell $j at time $(integrator.t)")
            terminate!(integrator, :InfDetected)
        end
    end
    update_values!(integrator.u, integrator.p, integrator.t)
end

function update_values!(U, params, t = 0)

    #update useful quantities relevant for potential, electron energy and fluid solve
    ncells = size(U, 2) - 2

    (;z_cell, fluids, fluid_ranges, index, scheme, source_term!, z_edge) = params
    (;BC_L, BC_R, B, ue, Tev, ∇ϕ, ϕ, pe, ne, νan, νc, μ) = params.cache
    OVS = params.OVS.energy.active

    mi = params.propellant.m

    # Apply fluid boundary conditions
    @views left_boundary_state!(BC_L, U[:, 2], params)
    @views right_boundary_state!(BC_R, U[:, end-1], params)

    @views U[:, 1] .= BC_L
    @views U[:, end] .= BC_R

    # Update electron quantities

    ne[1] = (1 - OVS) * max(params.config.min_number_density, electron_density(BC_L, params) / mi) + params.OVS.energy.ne(z_cell[1])
    ne[end] = (1 - OVS) * max(params.config.min_number_density, electron_density(BC_R, params) / mi) + params.OVS.energy.ne(z_cell[end])
    Tev[1] = (1 - OVS) * BC_L[index.nϵ]/ne[1] + params.OVS.energy.Tev(z_cell[1])
    Tev[end] = (1 - OVS) * BC_R[index.nϵ]/ne[end] + params.OVS.energy.Tev(z_cell[end])
    pe[1] = BC_L[index.nϵ]
    pe[end] = BC_R[index.nϵ]
    νan[1] = params.anom_model(BC_L, params, 1)
    νan[end] = params.anom_model(BC_R, params, length(z_cell))
    νc[1] = electron_collision_freq(Tev[1], BC_L[index.ρn]/mi , ne[1], mi)
    νc[end] = electron_collision_freq(Tev[end], BC_R[index.ρn]/mi , ne[end], mi)
    μ[1] = (1 - OVS)*electron_mobility(νan[1], νc[1], B[1])
    μ[end] = (1 - OVS)*electron_mobility(νan[end], νc[end], B[end])

    @inbounds for i in 2:(ncells + 1)
        z = z_cell[i]
        OVS_ne = OVS * (params.OVS.energy.ne(z))
        OVS_Tev = OVS * (params.OVS.energy.Tev(z))

        @views ne[i] = (1 - OVS) * max(params.config.min_number_density, electron_density(U[:, i], params) / mi) + OVS_ne
        Tev[i] = (1 - OVS) * U[index.nϵ, i]/ne[i] + OVS_Tev
        pe[i] = U[index.nϵ, i]
        @views params.cache.νan[i] = params.anom_model(U[:, i], params, i)
        params.cache.νc[i] = electron_collision_freq(params.cache.Tev[i], U[1, i]/mi , ne[i], mi)
        params.cache.μ[i] = (1 - params.OVS.energy.active)*electron_mobility(params.cache.νan[i], params.cache.νc[i], B[i]) #+ OVS*(params.OVS.energy.μ)
    end

    # update electrostatic potential and potential gradient on edges
    solve_potential_edge!(U, params)

    ∇ϕ[1] = forward_difference(ϕ[1], ϕ[2], ϕ[3], z_edge[1], z_edge[2], z_edge[3])
    ∇ϕ[end] = backward_difference(ϕ[end-2], ϕ[end-1], ϕ[end], z_edge[end-2], z_edge[end-1], z_edge[end])

    # Compute interior potential gradient and electron velocity and update source terms
    @inbounds for i in 2:(ncells + 1)
        # potential gradient
        ∇ϕ[i] = first_deriv_central_diff_pot(ϕ, params.z_cell, i)

        # electron velocity
        ue[i] = (1 - OVS) * electron_velocity(U, params, i) + OVS * (params.OVS.energy.ue)
    end

    ue[1] = ue[2] = ue[3]#electron_velocity(U, params, 1)
    ue[end] = ue[end-1] = ue[end-2] #electron_velocity(U, params, ncells+2)

    if params.implicit_energy > 0
        update_electron_energy_implicit!(U, params)
    else
        # Dirchlet BCs for electron energy
        apply_bc_electron!(U, params.BCs[3], :left, index)
        apply_bc_electron!(U, params.BCs[4], :right, index)
    end
end
