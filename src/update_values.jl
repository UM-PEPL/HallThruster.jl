
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
    (;BC_L, BC_R, B, ue, Tev, ∇ϕ, ϕ, pe, ne, νan, νc, μ, ∇pe) = params.cache
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

        @views ne[i] = (1 - OVS) * electron_density(U[:, i], params) / mi + OVS_ne
        Tev[i] = (1 - OVS) * U[index.nϵ, i]/ne[i] + OVS_Tev
        pe[i] = U[index.nϵ, i]
        @views params.cache.νan[i] = params.anom_model(U[:, i], params, i)
        params.cache.νc[i] = electron_collision_freq(params.cache.Tev[i], U[1, i]/mi , ne[i], mi)
        params.cache.μ[i] = (1 - OVS) *electron_mobility(params.cache.νan[i], params.cache.νc[i], B[i]) #+ OVS*(params.OVS.energy.μ)
    end

    # update electrostatic potential and potential gradient on edges
    solve_potential!(ϕ, U, params)

    # Compute potential gradient, pressure gradient, and electron velocity
    compute_gradients!(∇ϕ, ∇pe, ue, U, params)

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

    limiter(f) = sin(π*f)

    # Centered difference in interior cells
    @inbounds for j in 2:ncells-1
        cL, c0, cR = central_diff_coeffs(z_cell[j-1], z_cell[j], z_cell[j+1])
        for (f, g) in zip(functions, gradients)
            g[j] = cL * f[j-1] + c0 * f[j] + cR * f[j+1]
        end

        # Compute electron velocity (with limited gradients)
        f_ϕ = max(0.0, min(1.0, (ϕ[j] - ϕ[j-1])/(ϕ[j+1] - ϕ[j-1])))
        f_pe = max(0.0, min(1.0, (pe[j] - pe[j-1])/(pe[j+1] - pe[j-1])))
        ψ_ϕ = limiter(f_ϕ)
        ψ_pe = limiter(f_pe)

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