
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
    (;B, ue, Tev, ∇ϕ, ϕ, pe, ne, νan, νc, μ, ∇pe) = params.cache
    OVS = params.OVS.energy.active

    mi = params.propellant.m

    # Apply fluid boundary conditions
    @views left_boundary_state!(U[:, 1], U[:, 2], params)
    @views right_boundary_state!(U[:, end], U[:, end-1], params)

    # Update electron quantities
    @inbounds for i in 1:(ncells + 2)
        z = z_cell[i]
        OVS_ne = OVS * (params.OVS.energy.ne(z))
        OVS_Tev = OVS * (params.OVS.energy.Tev(z))

        @views ne[i] = (1 - OVS) * electron_density(U[:, i], params) / mi + OVS_ne
        Tev[i] = (1 - OVS) * U[index.nϵ, i]/ne[i] + OVS_Tev
        pe[i] = U[index.nϵ, i]
        @views params.cache.νan[i] = params.anom_model(U, params, i)
        params.cache.νc[i] = electron_collision_freq(params.cache.Tev[i], U[1, i]/mi , ne[i], mi)
        params.cache.μ[i] = (1 - OVS) * electron_mobility(params.cache.νan[i], params.cache.νc[i], B[i]) #+ OVS*(params.OVS.energy.μ)
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

    # Update electron energy if implicit, or if not then set electron boundary conditions for explicit solve
    if params.implicit_energy > 0
        update_electron_energy_implicit!(U, params)
    else
        # Dirchlet BCs for electron energy
        apply_bc_electron!(U, params.BCs[3], :left, index)
        apply_bc_electron!(U, params.BCs[4], :right, index)
    end
end

function compute_gradients!(∇ϕ, ∇pe, ue, U, params)
    (; ϕ, μ, ne, pe) = params.cache
    (;z_cell) = params

    ncells = length(z_cell)

    functions = (ϕ, pe)
    gradients = (∇ϕ, ∇pe)

    # Forward difference at left boundary
    cL, c0, cR = forward_diff_coeffs(z_cell[1], z_cell[2], z_cell[3])
    @inbounds for (f, g) in zip(functions, gradients)
        g[1] = cL * f[1] + c0 * f[2] + cR * f[3]
    end

    # Compute electron velocity
    ue[1] = μ[1] * (∇ϕ[1] - ∇pe[1]/ne[1])

    # Centered difference in interior cells
    @inbounds for j in 2:ncells-1
        cL, c0, cR = central_diff_coeffs(z_cell[j-1], z_cell[j], z_cell[j+1])
        for (f, g) in zip(functions, gradients)
            g[j] = cL * f[j-1] + c0 * f[j] + cR * f[j+1]
        end

        # Compute electron velocity
        ue[j] = μ[j] * (∇ϕ[j] - ∇pe[j]/ne[j])
    end

    # Backward difference at right boundary
    cL, c0, cR = backward_diff_coeffs(z_cell[end-2], z_cell[end-1], z_cell[end])
    @inbounds for (f, g) in zip(functions, gradients)
        g[end] = cL * f[end-2] + c0 * f[end-1] + cR * f[end]
    end

    # Compute electron velocity
    ue[end] = μ[end] * (∇ϕ[end] - ∇pe[end]/ne[end])

    return nothing
end