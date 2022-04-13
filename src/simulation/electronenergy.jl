function update_electron_energy!(U, params)
    (;Aϵ, bϵ, μ, ue, ne) = params.cache
    (;z_cell, dt, index) = params
    implicit = params.config.implicit_energy
    explicit = 1 - implicit
    ncells = size(U, 2)

    nϵ = @views U[index.nϵ, :]
    Aϵ.d[1] = 1.0
    Aϵ.du[1] = 0.0
    Aϵ.d[end] = 1.0
    Aϵ.dl[end] = 0.0

    bϵ[1] = 1.5 * params.Te_L * ne[1]
    bϵ[end] = 1.5 * params.Te_R * ne[end]

    @inbounds for i in 2:ncells-1
        Q = source_electron_energy(U, params, i)
        # User-provided source term
        Q += params.config.source_energy(U, params, i)

        zL = z_cell[i-1]
        z0 = z_cell[i]
        zR = z_cell[i+1]

        neL = ne[i-1]
        ne0 = ne[i]
        neR = ne[i+1]

        nϵL = nϵ[i-1]
        nϵ0 = nϵ[i]
        nϵR = nϵ[i+1]

        ϵL = nϵL / neL
        ϵ0 = nϵ0 / ne0
        ϵR = nϵR / neR

        ueL = ue[i-1]
        ue0 = ue[i]
        ueR = ue[i+1]

        μnϵL = μ[i-1] * nϵL
        μnϵ0 = μ[i] * nϵ0
        μnϵR = μ[i+1] * nϵR

        if params.config.LANDMARK
            # Use simplified thermal condutivity
            κL = μnϵL
            κ0 = μnϵ0
            κR = μnϵR
        else
            # Adjust thermal conductivity to be slightly more accurate
            κL = 24/25 * (1 / (1 + params.cache.νei[i-1] / √(2) / params.cache.νc[i-1])) * μnϵL
            κ0 = 24/25 * (1 / (1 + params.cache.νei[i]   / √(2) / params.cache.νc[i]))   * μnϵ0
            κR = 24/25 * (1 / (1 + params.cache.νei[i+1] / √(2) / params.cache.νc[i+1])) * μnϵR
        end

        # coefficients for centered three-point finite difference stencils
        d_cL, d_c0, d_cR = central_diff_coeffs(zL, z0, zR)
        d2_cL, d2_c0, d2_cR = second_deriv_coeffs(zL, z0, zR)

        # finite difference derivatives
        ∇nϵue = d_cL * ueL * nϵL + d_c0 * ue0 * nϵ0 + d_cR * ueR * nϵR
        ∇κ  = d_cL * κL  + d_c0 * κ0  + d_cR * κR
        ∇ϵ  = d_cL * ϵL  + d_c0 * ϵ0  + d_cR * ϵR
        ∇²ϵ = d2_cL * ϵL + d2_c0 * ϵ0 + d2_cR * ϵR

        # Explicit flux term
        F_explicit = 5/3 * ∇nϵue - 10/9 * (κ0 * ∇²ϵ + ∇κ * ∇ϵ)

        # Contribution to implicit part from μnϵ * d²ϵ/dz² term
        Aϵ.d[i]    = -10/9 * κ0 * d2_c0 / ne0
        Aϵ.dl[i-1] = -10/9 * κ0 * d2_cL / neL
        Aϵ.du[i]   = -10/9 * κ0 * d2_cR / neR

        # Contribution to implicit part from dμnϵ/dz * ∇ϵ term
        Aϵ.d[i]    -= 10/9 * ∇κ / ne0 * d_c0
        Aϵ.dl[i-1] -= 10/9 * ∇κ / neL * d_cL
        Aϵ.du[i]   -= 10/9 * ∇κ / neR * d_cR

        # Contribution to implicit part from advection term
        Aϵ.d[i]    += 5/3 * ue0 * d_c0
        Aϵ.dl[i-1] += 5/3 * ueL * d_cL
        Aϵ.du[i]   += 5/3 * ueR * d_cR

        # Contribution to implicit part from timestepping
        Aϵ.d[i]    = 1.0 + implicit * dt * Aϵ.d[i]
        Aϵ.dl[i-1] = implicit * dt * Aϵ.dl[i-1]
        Aϵ.du[i]   = implicit * dt * Aϵ.du[i]

        # Explicit right-hand-side
        bϵ[i] = nϵ[i] + dt * (Q - explicit * F_explicit)
    end

    # Solve equation system using Thomas algorithm
    tridiagonal_solve!(nϵ, Aϵ, bϵ)

    # Make sure Tev is positive, limit if below user-configured minumum electron temperature
    for i in 2:ncells-1
        if isnan(nϵ[i]) || isinf(nϵ[i]) || nϵ[i] / ne[i] < params.config.min_electron_temperature || nϵ[i] < params.config.min_electron_temperature
            nϵ[i] = params.config.min_electron_temperature * ne[i]
        end
    end
end