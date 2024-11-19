
function update_electron_energy!(params, dt)
    (; Δz_cell, Δz_edge, config, cache, Te_L, Te_R, ncells) = params
    (; A_energy, b_energy, nϵ, ue, ne, Tev, channel_area, dA_dz, κ, ni, niui) = cache
    implicit = params.config.implicit_energy
    explicit = 1 - implicit

    A_energy.d[1] = 1.0
    A_energy.du[1] = 0.0
    A_energy.d[end] = 1.0
    A_energy.dl[end] = 0.0

    if config.anode_boundary_condition == :dirichlet || ue[1] > 0
        b_energy[1] = 1.5 * Te_L * ne[1]
    else
        # Neumann BC for electron temperature
        b_energy[1] = 0
        A_energy.d[1] = 1.0 / ne[1]
        A_energy.du[1] = -1.0 / ne[2]
    end

    b_energy[end] = 1.5 * Te_R * ne[end]

    # Compute energy source terms
    Q = cache.cell_cache_1
    source_electron_energy!(Q, params)

    # Apply user-provided source terms
    if !isnothing(config.source_energy)
        @inbounds for i in 2:(ncells - 1)
            Q[i] += config.source_energy(params, i)
        end
    end

    @inbounds for i in 2:(ncells - 1)
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

        ΔzL = Δz_edge[left_edge(i)]
        ΔzR = Δz_edge[right_edge(i)]
        Δz = Δz_cell[i]

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
            if i == 2 && config.anode_boundary_condition == :sheath
                # left flux is sheath heat flux
                Te0 = 2 / 3 * nϵ0 / ne0

                # discharge current density
                jd = params.cache.Id[] / channel_area[1]

                # current densities at sheath edge
                ji_sheath_edge = e * sum(Z * niui[Z, 1] for Z in 1:(params.config.ncharge))
                je_sheath_edge = jd - ji_sheath_edge

                ne_sheath_edge = sum(Z * ni[Z, 1] for Z in 1:(params.config.ncharge))
                ue_sheath_edge = -je_sheath_edge / ne_sheath_edge / e

                FL_factor_L = 0.0
                FL_factor_C = 4 / 3 * ue_sheath_edge * (1 + params.cache.Vs[] / Te0)
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
        A_energy.d[i] = (FR_factor_C - FL_factor_C) / Δz
        A_energy.dl[i - 1] = (FR_factor_L - FL_factor_L) / Δz
        A_energy.du[i] = (FR_factor_R - FL_factor_R) / Δz

        # Contribution to implicit part from timestepping
        A_energy.d[i] = 1.0 + implicit * dt * A_energy.d[i]
        A_energy.dl[i - 1] = implicit * dt * A_energy.dl[i - 1]
        A_energy.du[i] = implicit * dt * A_energy.du[i]

        # Explicit flux
        F_explicit = (FR - FL) / Δz

        # Term to allow for changing area
        dlnA_dz = dA_dz[i] / channel_area[i]
        flux = 5 / 3 * nϵ0 * ue0

        # Explicit right-hand-side
        b_energy[i] = nϵ[i] + dt * (Q[i] - explicit * F_explicit)
        b_energy[i] -= dt * flux * dlnA_dz
    end

    # Solve equation system using Thomas algorithm
    tridiagonal_solve!(nϵ, A_energy, b_energy)

    # Make sure Tev is positive, limit if below user-configured minumum electron temperature
    @inbounds for i in 1:ncells
        if !isfinite(nϵ[i]) ||
           nϵ[i] / ne[i] < params.config.min_electron_temperature
            nϵ[i] = 3 / 2 * params.config.min_electron_temperature * ne[i]
        end

        # update Tev and pe
        Tev[i] = 2 / 3 * nϵ[i] / ne[i]
        if params.config.LANDMARK
            params.cache.pe[i] = nϵ[i]
        else
            params.cache.pe[i] = 2 / 3 * nϵ[i]
        end
    end
end
