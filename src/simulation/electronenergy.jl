const ELECTRON_CONDUCTIVITY_LOOKUP = let
    Z = [1.0, 2.0, 3.0, 4.0, 5.0]
    cond_coeff = [4.66, 4.0, 3.7, 3.6, 3.2]
    LinearInterpolation(Z, cond_coeff)
end

function update_electron_energy!(U, params)
    (;z_cell, z_edge, dt, index, config, cache, Te_L, Te_R) = params
    (;Aϵ, bϵ, μ, ue, ne, Tev) = cache
    implicit = params.config.implicit_energy
    explicit = 1 - implicit
    ncells = size(U, 2)

    nϵ = @views U[index.nϵ, :]
    Aϵ.d[1] = 1.0
    Aϵ.du[1] = 0.0
    Aϵ.d[end] = 1.0
    Aϵ.dl[end] = 0.0

    if config.anode_boundary_condition == :dirichlet || ue[1] > 0
        bϵ[1] = 1.5 * Te_L * ne[1]
    else
        # Neumann BC for internal energy
        bϵ[1] = 0
        Aϵ.d[1] = 1.0
        Aϵ.du[1] = -1.0
    end

    bϵ[end] = 1.5 * Te_R * ne[end]

    Δt = dt

    # Needed to compute excitation and ionization frequencies in first and last cells,
    # Need a better solution, because the signature of this function doesn't make it clear
    # That params.cache.νex and params.cache.νei are being modified
    _ = inelastic_losses!(U, params, 1)
    _ = inelastic_losses!(U, params, ncells)

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

        μnϵL = μ[i-1] * nϵL
        μnϵ0 = μ[i] * nϵ0
        μnϵR = μ[i+1] * nϵR

        if config.LANDMARK
            # Use simplified thermal condutivity
            κL = 10/9 * μnϵL
            κ0 = 10/9 * μnϵ0
            κR = 10/9 * μnϵR

        else
            #get adjusted coeffient for higher charge states
            κ_charge = ELECTRON_CONDUCTIVITY_LOOKUP(params.cache.Z_eff[i])
            correction_factor = κ_charge/4.7
            # Adjust thermal conductivity to be slightly more accurate
            κL = 10/9 * 24/25 * (1 / (1 + params.cache.νei[i-1] / √(2) / params.cache.νc[i-1])) * μnϵL * correction_factor
            κ0 = 10/9 * 24/25 * (1 / (1 + params.cache.νei[i]   / √(2) / params.cache.νc[i]))   * μnϵ0 * correction_factor
            κR = 10/9 * 24/25 * (1 / (1 + params.cache.νei[i+1] / √(2) / params.cache.νc[i+1])) * μnϵR * correction_factor
        end

        # Weighted average of the electron velocities in the three stencil cells
        ue_avg = 0.25 * (ueL + 2 * ue0 + ueR)

        ΔzL = z_cell[i] - z_cell[i-1]
        ΔzR = z_cell[i+1] - z_cell[i]
        Δz = z_edge[right_edge(i)] - z_edge[left_edge(i)]

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

                uth = -0.25 * sqrt(8 * e * Te0 / π / me) * exp(-params.cache.Vs[] / Te0)

                FL_factor_L = 0.0
                FL_factor_C = 4/3 * uth
                FL_factor_R = 0.0

                Q += ne0 * uth * params.cache.Vs[] / Δz
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
        Aϵ.d[i]    = 1.0 + implicit * Δt * Aϵ.d[i]
        Aϵ.dl[i-1] = implicit * Δt * Aϵ.dl[i-1]
        Aϵ.du[i]   = implicit * Δt * Aϵ.du[i]

        # Explicit flux
        F_explicit = (FR - FL) / Δz

        # Explicit right-hand-side
        bϵ[i] = nϵ[i] + Δt * (Q - explicit * F_explicit)
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
