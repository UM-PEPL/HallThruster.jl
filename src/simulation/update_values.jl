
function update_values!(integrator)
    (nvars, ncells) = size(integrator.u)

    nandetected = false
    infdetected = false

    @inbounds for j in 1:ncells, i in 1:nvars
        if isnan(integrator.u[i, j])
            println("NaN detected in variable $i in cell $j at time $(integrator.t)")
            nandetected = true
            terminate!(integrator, :NaNDetected)
            break
        elseif isinf(integrator.u[i, j])
            println("Inf detected in variable $i in cell $j at time $(integrator.t)")
            infdetected = true
            terminate!(integrator, :InfDetected)
            break
        end
    end

    if !nandetected && !infdetected
        update_values!(integrator.u, integrator.p, integrator.t)
    end
end

#update useful quantities relevant for potential, electron energy and fluid solve
function update_values!(U, params, t = 0)
    (;z_cell, index, num_subiterations) = params
    (;B, ue, Tev, ∇ϕ, ϕ, pe, ne, μ, ∇pe, νan, νc, νen, νei, νw, Z_eff, νe) = params.cache

    mi = params.config.propellant.m

    # Update the current iteration
    params.iteration[1] += 1

    # Apply fluid boundary conditions
    @views left_boundary_state!(U[:, 1], U, params)
    @views right_boundary_state!(U[:, end], U, params)

    ncells = size(U, 2) - 2

    for i in 1:num_subiterations
        # Update electron quantities
        @inbounds for i in 1:(ncells + 2)
            z = z_cell[i]

            ne[i] = max(params.config.min_number_density, electron_density(U, params, i))
            Tev[i] = 2/3 * max(params.config.min_electron_temperature, U[index.nϵ, i]/ne[i])
            pe[i] = if params.config.LANDMARK
                3/2 * ne[i] * Tev[i]
            else
                ne[i] * Tev[i]
            end #U[index.nϵ, i]
            νen[i] = freq_electron_neutral(U, params, i)
            νei[i] = freq_electron_ion(U, params, i)
            νw[i] = freq_electron_wall(U, params, i)
            νan[i] = freq_electron_anom(U, params, i)
            νc[i] = νen[i] + νei[i]
            νe[i] = νc[i] + νan[i] + νw[i]
            μ[i] = electron_mobility(νe[i], B[i])
            Z_eff[i] = compute_Z_eff(U, params, i)
        end

        # update electrostatic potential and potential gradient on edges
        solve_potential_edge!(U, params)
        # Compute potential gradient, pressure gradient, and electron velocity
        compute_gradients!(∇ϕ, ∇pe, ue, U, params)

        # Fix electron velocity on left and right cells
        ueL = ue[3]
        ue[1] = ue[2] = ueL
        ueR = ue[end-2]
        ue[end] = ue[end-1] = ueR

        update_electron_energy!(U, params)
    end

end

function compute_gradients!(∇ϕ, ∇pe, ue, U, params)
    (; ϕ, μ, ne, pe, ϕ_cell) = params.cache
    (;z_cell, z_edge) = params

    ncells = length(z_cell)

    # Interpolate potential to cells
    ϕ_cell[1] = ϕ[1]
    ϕ_cell[end] = ϕ[end]
    @turbo for i in 2:ncells-1
        ϕ_cell[i] = lerp(z_cell[i], z_edge[i-1], z_edge[i], ϕ[i-1], ϕ[i])
    end

    # Potential gradient (centered)
    ∇ϕ[1] = forward_difference(ϕ[1], ϕ[2], ϕ[3], z_edge[1], z_edge[2], z_edge[3])
    #∇ϕ[1] = (ϕ[2] - ϕ[1]) / (z_edge[2] - z_edge[1])
    # Pressure gradient (forward)
    ∇pe[1] = forward_difference(pe[1], pe[2], pe[3], z_cell[1], z_cell[2], z_cell[3])
    # Compute electron velocity
    ue[1] = μ[1] * (∇ϕ[1] - ∇pe[1]/ne[1])

    # Centered difference in interior cells
    @inbounds for j in 2:ncells-1

        # Compute potential gradient
        ∇ϕ[j] = (ϕ[j] - ϕ[j-1]) / (z_edge[j] - z_edge[j-1])

        # Compute pressure gradient
        ∇pe[j] = central_difference(pe[j-1], pe[j], pe[j+1], z_cell[j-1], z_cell[j], z_cell[j+1])
        # Compute electron velocity
        ue[j] = μ[j] * (∇ϕ[j] - ∇pe[j]/ne[j])
    end

    # Potential gradient (centered)
    ∇ϕ[end] = backward_difference(ϕ[end-2], ϕ[end-1], ϕ[end], z_edge[end-2], z_edge[end-1], z_edge[end])
    #∇ϕ[end] = (ϕ[end] - ϕ[end-1]) / (z_edge[end] - z_edge[end-1])
    # pressure gradient (backward)
    ∇pe[end] = backward_difference(pe[end-2], pe[end-1], pe[end], z_cell[end-2], z_cell[end-1], z_cell[end])
    # Compute electron velocity
    ue[end] = μ[end] * (∇ϕ[end] - ∇pe[end]/ne[end])

    return nothing
end