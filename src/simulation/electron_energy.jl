function update_electron_energy!(params, wall_loss_model, source_energy, dt)
    (; Te_L, Te_R, grid, cache, min_Te, implicit_energy, anode_bc, landmark) = params
    (; Aϵ, bϵ, nϵ, ue, ne, Tev, pe) = cache

    Q = cache.cell_cache_1

    # Compute energy source terms
    source_electron_energy!(Q, params, wall_loss_model)

    # User-provided source term
    @inbounds for i in 2:(grid.num_cells - 1)
        Q[i] += source_energy(params, i)
    end

    energy_boundary_conditions!(Aϵ, bϵ, Te_L, Te_R, ne, ue, anode_bc)

    # Fill matrix and RHS
    setup_energy_system!(Aϵ, bϵ, grid, cache, anode_bc, implicit_energy, dt)

    # Solve equation system using the Thomas algorithm
    tridiagonal_solve!(nϵ, Aϵ, bϵ)

    # Make sure Tev is positive, limit if below minumum electron temperature
    limit_energy!(nϵ, ne, min_Te)
    update_temperature!(Tev, nϵ, ne, min_Te)
    update_pressure!(pe, nϵ, landmark)
    return
end

function setup_energy_system!(Aϵ, bϵ, grid, cache, anode_bc, implicit, dt)
    (; ne, ue, κ, ji, channel_area, dA_dz, nϵ) = cache

    explicit = 1.0 - implicit
    Q = cache.cell_cache_1

    @inbounds for i in 2:(grid.num_cells - 1)
        neL = ne[i - 1]
        ne0 = ne[i]
        neR = ne[i + 1]

        nϵL = nϵ[i - 1]
        nϵ0 = nϵ[i]
        nϵR = nϵ[i + 1]

        ueL = ue[i - 1]
        ue0 = ue[i]
        ueR = ue[i + 1]

        #pull thermal conductivity values
        κL = κ[i - 1]
        κ0 = κ[i]
        κR = κ[i + 1]

        ΔzL = grid.dz_edge[left_edge(i)]
        ΔzR = grid.dz_edge[right_edge(i)]
        Δz = grid.dz_cell[i]

        # Weighted average of the electron velocities in the three stencil cells
        ue_avg = 0.25 * (ΔzL * (ueL + ue0) + ΔzR * (ue0 + ueR)) / Δz

        # Upwind differences
        if ue_avg > 0
            FR_factor_L = 0.0
            FR_factor_C = 5 / 3 * ue0 + κ0 / ΔzR / ne0
            FR_factor_R = -κ0 / ΔzR / neR

            FL_factor_L = 5 / 3 * ueL + κL / ΔzL / neL
            FL_factor_C = -κL / ΔzL / ne0
            FL_factor_R = 0.0
        else
            if i == 2 && anode_bc == :sheath
                # left flux is sheath heat flux
                Te0 = 2 / 3 * nϵ0 / ne0

                # discharge current density
                jd = cache.Id[] / channel_area[1]

                # current densities at sheath edge
                ji_sheath_edge = ji[1]
                je_sheath_edge = jd - ji_sheath_edge

                ne_sheath_edge = ne[1]
                ue_sheath_edge = -je_sheath_edge / ne_sheath_edge / e

                FL_factor_L = 0.0
                FL_factor_C = 4 / 3 * ue_sheath_edge * (1 + cache.Vs[] / Te0)
                FL_factor_R = 0.0
            elseif i == 2
                # central differences at left boundary for compatibility with dirichlet BC
                FL_factor_L = 5 / 3 * ueL + κL / ΔzL / neL
                FL_factor_C = -κL / ΔzL / ne0
                FL_factor_R = 0.0
            else
                # Upwind differences
                FL_factor_L = κ0 / ΔzL / neL
                FL_factor_C = 5 / 3 * ue0 - κ0 / ΔzL / ne0
                FL_factor_R = 0.0
            end

            FR_factor_L = 0.0
            FR_factor_C = κR / ΔzR / ne0
            FR_factor_R = 5 / 3 * ueR - κR / ΔzR / neR
        end

        # Fluxes at left and right boundary
        FL = FL_factor_L * nϵL + FL_factor_C * nϵ0 + FL_factor_R * nϵR
        FR = FR_factor_L * nϵL + FR_factor_C * nϵ0 + FR_factor_R * nϵR

        # Contribution to implicit part from fluxes
        Aϵ.d[i] = (FR_factor_C - FL_factor_C) / Δz
        Aϵ.dl[i - 1] = (FR_factor_L - FL_factor_L) / Δz
        Aϵ.du[i] = (FR_factor_R - FL_factor_R) / Δz

        # Contribution to implicit part from timestepping
        Aϵ.d[i] = 1.0 + implicit * dt * Aϵ.d[i]
        Aϵ.dl[i - 1] = implicit * dt * Aϵ.dl[i - 1]
        Aϵ.du[i] = implicit * dt * Aϵ.du[i]

        # Explicit flux
        F_explicit = (FR - FL) / Δz

        # Term to allow for changing area
        dlnA_dz = dA_dz[i] / channel_area[i]
        flux = 5 / 3 * nϵ0 * ue0

        # Explicit right-hand-side
        bϵ[i] = nϵ[i] + dt * (Q[i] - explicit * F_explicit)
        bϵ[i] -= dt * flux * dlnA_dz
    end
    return
end

function energy_boundary_conditions!(Aϵ, bϵ, Te_L, Te_R, ne, ue, anode_bc)
    Aϵ.d[1] = 1.0
    Aϵ.du[1] = 0.0
    Aϵ.d[end] = 1.0
    Aϵ.dl[end] = 0.0

    if anode_bc == :dirichlet || ue[1] > 0
        bϵ[1] = 1.5 * Te_L * ne[1]
    else
        # Neumann BC for electron temperature
        bϵ[1] = 0
        Aϵ.d[1] = 1.0 / ne[1]
        Aϵ.du[1] = -1.0 / ne[2]
    end

    bϵ[end] = 1.5 * Te_R * ne[end]
    return
end

function limit_energy!(nϵ, ne, min_Te)
    @inbounds for i in eachindex(nϵ)
        if !isfinite(nϵ[i]) || nϵ[i] / ne[i] < 1.5 * min_Te
            nϵ[i] = 1.5 * min_Te * ne[i]
        end
    end
    return
end

function update_temperature!(Tev, nϵ, ne, min_Te)
    # Update plasma quantities based on new density
    @inbounds for i in eachindex(Tev)
        # Compute new electron temperature
        Tev[i] = 2.0 / 3.0 * max(min_Te, nϵ[i] / ne[i])
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

#===============================================================================
Electron energy source terms
===============================================================================#

function excitation_losses!(Q, cache, landmark, grid, reactions, reactant_indices, fluids)
    (; νex, ϵ, ne, K) = cache
    ncells = length(grid.cell_centers)

    @. νex = 0.0
    for (ind, rxn) in zip(reactant_indices, reactions)
        dens = fluids[ind].density
        inv_m = 1 / fluids[ind].species.element.m
        for i in 2:(ncells - 1)
            r = rate_coeff(rxn, ϵ[i])
            ndot = reaction_rate(r, ne[i], dens[i] * inv_m)
            νex[i] += ndot / ne[i]
            Q[i] += ndot * (rxn.energy - !landmark * K[i])
        end
    end

    return nothing
end

function ohmic_heating!(Q, cache, landmark)
    (; ne, ue, ∇ϕ, K, νe, ue, ∇pe) = cache
    # Compute ohmic heating term, which is the rate at which energy is transferred out of the electron
    # drift (kinetic energy) into thermal energy
    if (landmark)
        # Neglect kinetic energy, so the rate of energy transfer into thermal energy is equivalent to
        # the total input power into the electrons (j⃗ ⋅ E⃗ = -mₑnₑ|uₑ|²νₑ)
        @. Q = ne * ue * ∇ϕ
    else
        # Do not neglect kinetic energy, so ohmic heating term is mₑnₑ|uₑ|²νₑ + ue ∇pe = 2nₑKνₑ + ue ∇pe
        # where K is the electron bulk kinetic energy, 1/2 * mₑ|uₑ|²
        @. Q = 2 * ne * K * νe + ue * ∇pe
    end
    return nothing
end

function source_electron_energy!(Q, params, wall_loss_model)
    (; cache, landmark, grid, excitation_reactions) = params
    (; ne, ohmic_heating, wall_losses, inelastic_losses) = cache

    # compute ohmic heating
    ohmic_heating!(ohmic_heating, cache, landmark)

    # add excitation losses to total inelastic losses
    excitation_losses!(
        inelastic_losses, cache, landmark, grid,
        excitation_reactions, params.excitation_reactant_indices,
        params.fluid_array
    )

    # compute wall losses
    wall_power_loss!(wall_losses, wall_loss_model, params)

    # Compute net energy source, i.e heating minus losses
    @. Q = ohmic_heating - ne * wall_losses - inelastic_losses

    return Q
end