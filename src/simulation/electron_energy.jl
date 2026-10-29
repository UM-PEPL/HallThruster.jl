function update_electron_energy!(params, wall_loss_model, user_source_energy!, dt)
    (; Te_L, Te_R, grid, cache, implicit_energy, anode_bc, landmark) = params
    (; Aϵ, bϵ, nϵ, ue, ne, Tev, pe, ∇pe, user_energy_source) = cache

    Q = cache.cell_cache_1

    # User-provided source term
    user_source_energy!(user_energy_source, params)

    # Compute other energy source terms
    source_electron_energy!(Q, params, wall_loss_model)

    # Set up energy boundary conditions
    energy_boundary_conditions!(Aϵ, bϵ, Te_L, Te_R, ne, ue, anode_bc)

    # Fill matrix and RHS
    setup_energy_system!(Aϵ, bϵ, grid, cache, anode_bc, implicit_energy, dt)

    # Solve equation system using the Thomas algorithm
    tridiagonal_solve!(nϵ, Aϵ, bϵ)

    # Limit energy and update temperature and pressure in one state traversal.
    update_energy_state!(Tev, nϵ, ne, pe, params.min_Te, landmark)
    update_pressure_gradient!(∇pe, pe, grid.cell_centers)

    return
end

# One pair of coefficients per face is shared by both adjacent cells. κ is
# defined against temperature; the unknown nϵ/ne is mean energy = 3/2 Tev.
function energy_flux_coefficients!(grid, cache, anode_bc)
    (; ne, ue, κ, nϵ, energy_flux_left, energy_flux_right) = cache
    @inbounds for face in eachindex(grid.edges)
        left, right = face, face + 1
        distance = grid.dz_edge[face]
        weight = (grid.edges[face] - grid.cell_centers[left]) / distance
        velocity = lerp(ue[left], ue[right], weight)
        conductivity = (2 / 3) * lerp(κ[left], κ[right], weight)
        energy_flux_left[face] = (5 / 3) * max(velocity, 0.0) +
            conductivity / (distance * ne[left])
        energy_flux_right[face] = (5 / 3) * min(velocity, 0.0) -
            conductivity / (distance * ne[right])
    end

    if anode_bc == :sheath && ue[2] <= 0
        # Maxwellian electrons carry 2 Tev + Vs out of the sheath edge.
        temperature = (2 / 3) * nϵ[2] / ne[2]
        energy_flux_left[1] = 0.0
        energy_flux_right[1] = (4 / 3) * ue[2] *
            (1 + cache.Vs[] / (2 * temperature))
    end
    return nothing
end

@inline function energy_face_area(grid, area, face)
    weight = (grid.edges[face] - grid.cell_centers[face]) / grid.dz_edge[face]
    return lerp(area[face], area[face + 1], weight)
end

function setup_energy_system!(Aϵ, bϵ, grid, cache, anode_bc, implicit, dt)
    (; channel_area, nϵ, energy_flux_left, energy_flux_right) = cache
    energy_flux_coefficients!(grid, cache, anode_bc)
    Q = cache.cell_cache_1

    @inbounds for i in interior_cells(grid.cell_centers)
        left, right = left_edge(i), right_edge(i)
        volume = channel_area[i] * grid.dz_cell[i]
        left_scale = energy_face_area(grid, channel_area, left) / volume
        right_scale = energy_face_area(grid, channel_area, right) / volume

        # A^-1 ∂z(A F) for the entire flux, including heat conduction. Sharing
        # face areas and fluxes makes internal transfers cancel exactly.
        lower = -left_scale * energy_flux_left[left]
        diagonal = right_scale * energy_flux_left[right] -
            left_scale * energy_flux_right[left]
        upper = right_scale * energy_flux_right[right]
        divergence = lower * nϵ[i - 1] + diagonal * nϵ[i] + upper * nϵ[i + 1]

        Aϵ.dl[i - 1] = implicit * dt * lower
        Aϵ.d[i] = 1.0 + implicit * dt * diagonal
        Aϵ.du[i] = implicit * dt * upper
        bϵ[i] = nϵ[i] + dt * (Q[i] - (1 - implicit) * divergence)
    end
    return nothing
end

function energy_boundary_conditions!(Aϵ, bϵ, Te_L, Te_R, ne, ue, anode_bc)
    if anode_bc == :dirichlet || ue[2] > 0
        # 0.5 (nϵ[1]/ne[1] + nϵ[2]/ne[2]) = 3/2 Te_L
        Aϵ.d[1] = 0.5 / ne[1]
        Aϵ.du[1] = 0.5 / ne[2]
        bϵ[1] = 1.5 * Te_L
    else
        # Neumann BC for electron temperature
        Aϵ.d[1] = 1.0 / ne[1]
        Aϵ.du[1] = -1.0 / ne[2]
        bϵ[1] = 0
    end

    # 0.5 (nϵ[end-1]/ne[end-1] + nϵ[end]/ne[end]) = 3/2 Te_R
    Aϵ.d[end] = 0.5 / ne[end]
    Aϵ.dl[end] = 0.5 / ne[end - 1]
    bϵ[end] = 1.5 * Te_R
    return
end

function update_energy_state!(Tev, nϵ, ne, pe, min_Te, landmark)
    all(isfinite, nϵ) || throw(DomainError(nϵ, "Nonfinite electron energy after the energy solve."))
    pe_factor = landmark ? 1.0 : 2.0 / 3.0
    @inbounds @simd for i in eachindex(Tev)
        energy = nϵ[i]
        temperature = max(min_Te, energy / ne[i] / 1.5)
        energy = 1.5 * ne[i] * temperature
        Tev[i] = temperature
        nϵ[i] = energy
        pe[i] = pe_factor * energy
    end
    return
end

function update_pressure!(pe, nϵ, landmark)
    if landmark
        pe_factor = 1.0
    else
        pe_factor = 2.0 / 3.0
    end
    @inbounds for i in eachindex(pe)
        pe[i] = pe_factor * nϵ[i]
    end
    return
end

function update_pressure_gradient!(∇pe, pe, z)
    # Pressure gradient (forward)
    ∇pe[1] = forward_difference(pe[1], pe[2], pe[3], z[1], z[2], z[3])

    # Centered difference in interior cells
    @inbounds for j in interior_cells(pe)
        # Compute pressure gradient
        ∇pe[j] = central_difference(pe[j - 1], pe[j], pe[j + 1], z[j - 1], z[j], z[j + 1])
    end

    # pressure gradient (backward)
    ∇pe[end] = backward_difference(pe[end - 2], pe[end - 1], pe[end], z[end - 2], z[end - 1], z[end])

    return nothing
end

#===============================================================================
Electron energy source terms
===============================================================================#

function excitation_losses!(Q, cache, reactions, reactant_indices, fluids)
    (; ϵ, ne) = cache

    for (ind, rxn) in zip(reactant_indices, reactions)
        dens = fluids[ind].density
        inv_m = 1 / fluids[ind].species.element.m
        @inbounds for i in interior_cells(Q)
            r = rate_coeff(rxn, ϵ[i])
            ndot = reaction_rate(r, ne[i], dens[i] * inv_m)
            Q[i] += ndot * rxn.energy
        end
    end

    return nothing
end

function add_lumped_excitation_frequency!(
        νex, cache, reactions, reactant_indices, fluids,
    )
    (; ϵ) = cache

    for (ind, rxn) in zip(reactant_indices, reactions)
        dens = fluids[ind].density
        inv_m = 1 / fluids[ind].species.element.m
        @inbounds for i in interior_cells(νex)
            r = rate_coeff(rxn, ϵ[i])
            excitation_frequency = r * dens[i] * inv_m
            νex[i] += excitation_frequency
        end
    end

    return nothing
end

function ohmic_heating!(Q, cache, landmark)
    (; ne, ue, ∇ϕ, K, νe, ∇pe) = cache
    # Compute ohmic heating term, which is the rate at which energy is transferred out of the electron
    # drift (kinetic energy) into thermal energy
    if (landmark)
        # Neglect kinetic energy, so the rate of energy transfer into thermal energy is equivalent to
        # the total input power into the electrons (j⃗ ⋅ E⃗ = -mₑnₑ|uₑ|²νₑ)
        @inbounds for i in 2:(length(Q) - 1)
            Q[i] = ne[i] * ue[i] * ∇ϕ[i]
        end
    else
        # Do not neglect kinetic energy, so ohmic heating term is mₑnₑ|uₑ|²νₑ + ue ∇pe = 2nₑKνₑ + ue ∇pe
        # where K is the electron bulk kinetic energy, 1/2 * mₑ|uₑ|²
        @inbounds for i in 2:(length(Q) - 1)
            Q[i] = 2 * ne[i] * K[i] * νe[i] + ue[i] * ∇pe[i]
        end
    end
    return nothing
end

function source_electron_energy!(Q, params, wall_loss_model)
    (; cache, landmark, excitation_reactions) = params
    (; ne, ohmic_heating, wall_losses, inelastic_losses, user_energy_source) = cache

    # compute ohmic heating
    ohmic_heating!(ohmic_heating, cache, landmark)

    # add excitation losses to total inelastic losses
    excitation_losses!(
        inelastic_losses, cache, excitation_reactions,
        params.excitation_reactant_indices, params.fluid_array,
    )

    # compute wall losses
    wall_power_loss!(wall_losses, wall_loss_model, params)

    # Compute net energy source, i.e heating minus losses
    @inbounds for i in 2:(params.grid.num_cells - 1)
        Q[i] = user_energy_source[i] + ohmic_heating[i] - ne[i] * wall_losses[i] - inelastic_losses[i]
    end

    return Q
end
