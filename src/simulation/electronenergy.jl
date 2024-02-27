

function update_electron_energy!(U, params, dt)
    (;Δz_cell, Δz_edge, index, config, cache, Te_L, Te_R) = params
    (;Aϵ, bϵ, nϵ, ue, ne, Tev, channel_area, dA_dz, κ) = cache
    implicit = params.config.implicit_energy
    explicit = 1 - implicit
    ncells = size(U, 2)
    mi = params.config.propellant.m

    Aϵ.d[1] = 1.0
    Aϵ.du[1] = 0.0
    Aϵ.d[end] = 1.0
    Aϵ.dl[end] = 0.0

    if config.anode_boundary_condition == :dirichlet || ue[1] > 0
        bϵ[1] = 1.5 * Te_L * ne[1]
    else
        # Neumann BC for electron temperature
        bϵ[1] = 0
        Aϵ.d[1] = 1.0 / ne[1]
        Aϵ.du[1] = -1.0 / ne[2]
    end

    bϵ[end] = 1.5 * Te_R * ne[end]

    @inbounds for i in 2:ncells-1
        Q = source_electron_energy(U, params, i)
        # User-provided source term
        Q += config.source_energy(U, params, i)

        neL = ne[i-1]
        ne0 = ne[i]
        neR = ne[i+1]

        nϵL = nϵ[i-1]
        nϵ0 = nϵ[i]
        nϵR = nϵ[i+1]

        ueL = ue[i-1]
        ue0 = ue[i]
        ueR = ue[i+1]

        #pull thermal conductivity values
        κL = κ[i-1]
        κ0 = κ[i]
        κR = κ[i+1]

        ΔzL = Δz_edge[left_edge(i)]
        ΔzR = Δz_edge[right_edge(i)]
        Δz = Δz_cell[i]

        # Weighted average of the electron velocities in the three stencil cells
        ue_avg = 0.25 * (ΔzL * (ueL + ue0) + ΔzR * (ue0 + ueR)) / Δz

        # Upwind differences
        if ue_avg > 0
            FR_factor_L = 0.0
            FR_factor_C = 5/3 * ue0 + κ0 / ΔzR / ne0
            FR_factor_R = - κ0 / ΔzR / neR

            FL_factor_L = 5/3 * ueL + κL / ΔzL / neL
            FL_factor_C = -κL / ΔzL / ne0
            FL_factor_R = 0.0
        else
            if i == 2 && config.anode_boundary_condition == :sheath
                # left flux is sheath heat flux
                Te0 = 2/3 * nϵ0 / ne0

                # discharge current density
                jd = params.cache.Id[] / channel_area[1]

                # current densities at sheath edge
                ji_sheath_edge = e * sum(Z * U[index.ρiui[Z], 1] for Z in 1:params.config.ncharge) / mi
                je_sheath_edge = jd - ji_sheath_edge

                ne_sheath_edge = sum(Z * U[index.ρi[Z], 1] for Z in 1:params.config.ncharge) / mi
                ue_sheath_edge = -je_sheath_edge / ne_sheath_edge / e

                FL_factor_L = 0.0
                FL_factor_C = 4/3 * ue_sheath_edge * (1 + params.cache.Vs[] / Te0)
                FL_factor_R = 0.0
            elseif i == 2
                # central differences at left boundary for compatibility with dirichlet BC
                FL_factor_L = 5/3 * ueL + κL / ΔzL / neL
                FL_factor_C = -κL / ΔzL / ne0
                FL_factor_R = 0.0
            else
                # Upwind differences
                FL_factor_L = κ0 / ΔzL / neL
                FL_factor_C = 5/3 * ue0 - κ0 / ΔzL / ne0
                FL_factor_R = 0.0
            end

            FR_factor_L = 0.0
            FR_factor_C = κR / ΔzR / ne0
            FR_factor_R = 5/3 * ueR - κR / ΔzR / neR
        end

        # Fluxes at left and right boundary
        FL = FL_factor_L * nϵL + FL_factor_C * nϵ0 + FL_factor_R * nϵR
        FR = FR_factor_L * nϵL + FR_factor_C * nϵ0 + FR_factor_R * nϵR

        # Contribution to implicit part from fluxes
        Aϵ.d[i] = (FR_factor_C - FL_factor_C) / Δz
        Aϵ.dl[i-1] = (FR_factor_L - FL_factor_L) / Δz
        Aϵ.du[i] = (FR_factor_R - FL_factor_R) / Δz

        # Contribution to implicit part from timestepping
        Aϵ.d[i]    = 1.0 + implicit * dt * Aϵ.d[i]
        Aϵ.dl[i-1] = implicit * dt * Aϵ.dl[i-1]
        Aϵ.du[i]   = implicit * dt * Aϵ.du[i]

        # Explicit flux
        F_explicit = (FR - FL) / Δz

        # Term to allow for changing area
        dlnA_dz = dA_dz[i] / channel_area[i]
        flux = 5/3 * nϵ0 * ue0

        # Explicit right-hand-side
        bϵ[i] = nϵ[i] + dt * (Q - explicit * F_explicit)
        bϵ[i] -= dt * flux * dlnA_dz
    end

    # Solve equation system using Thomas algorithm
    tridiagonal_solve!(nϵ, Aϵ, bϵ)

    # Make sure Tev is positive, limit if below user-configured minumum electron temperature
    @inbounds for i in 1:ncells

        if isnan(nϵ[i]) || isinf(nϵ[i]) || nϵ[i] / ne[i] < params.config.min_electron_temperature
            nϵ[i] = 3/2 *  params.config.min_electron_temperature * ne[i]
        end

        # update Tev and pe
        Tev[i] = 2/3 * nϵ[i] / ne[i]
        if params.config.LANDMARK
            params.cache.pe[i] = nϵ[i]
        else
            params.cache.pe[i] = 2/3 * nϵ[i]
        end
    end
end
