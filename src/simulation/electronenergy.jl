function ElectronCondLookup()
    Z = [1.0, 2.0, 3.0, 4.0, 5.0]
    cond_coeff = [4.66, 4.0, 3.7, 3.6, 3.2]
    coeff = LinearInterpolation(Z, cond_coeff)
    return coeff
end

function update_electron_energy!(U, params)
    (;Aϵ, bϵ, μ, ue, ne, Tev) = params.cache
    (;z_cell, dt, index) = params
    implicit = params.config.implicit_energy
    explicit = 1 - implicit
    ncells = size(U, 2)

    nϵ = @views U[index.nϵ, :]
    Aϵ.d[1] = 1.0
    Aϵ.du[1] = 0.0
    Aϵ.d[end] = 1.0
    Aϵ.dl[end] = 0.0

    #if params.config.LANDMARK
        bϵ[1] = 1.5 * params.Te_L * ne[1]
    #=else
        # Neumann BC for internal energy
        bϵ[1] = 0
        Aϵ.d[1] = 1.0
        Aϵ.du[1] = -1.0
    end=#

    bϵ[end] = 1.5 * params.Te_R * ne[end]

    Δt = dt

    # Needed to compute excitation and ionization frequencies in first and last cells,
    # Need a better solution, because the signature of this function doesn't make it clear
    # That params.cache.νex and params.cache.νei are being modified
    _ = source_electron_energy(U, params, 1)
    _ = source_electron_energy(U, params, ncells)

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

            flux_factor = 5/3
        else
            #get adjusted coeffient for higher charge states
            κ_charge = params.config.electron_cond_lookup(params.cache.Z_eff[i])
            correction_factor = κ_charge/4.7
            # Adjust thermal conductivity to be slightly more accurate
            κL = 24/25 * (1 / (1 + params.cache.νei[i-1] / √(2) / params.cache.νc[i-1])) * μnϵL * correction_factor
            κ0 = 24/25 * (1 / (1 + params.cache.νei[i]   / √(2) / params.cache.νc[i]))   * μnϵ0 * correction_factor
            κR = 24/25 * (1 / (1 + params.cache.νei[i+1] / √(2) / params.cache.νc[i+1])) * μnϵR * correction_factor
            flux_factor = 1.0
        end

        # coefficients for centered three-point finite difference stencils
        d_cL, d_c0, d_cR =
            if params.electron_energy_order == 2 || i == 2
                central_diff_coeffs(zL, z0, zR)
            else
                if ue[i] ≤ 0
                    downwind_diff_coeffs(zL, z0, zR)
                else
                    upwind_diff_coeffs(zL, z0, zR)
                end
            end
        d2_cL, d2_c0, d2_cR = second_deriv_coeffs(zL, z0, zR)

        Vs = 0.0

        #if params.config.LANDMARK || i > 2
            # finite difference derivatives
            ∇nϵue = d_cL * ueL * nϵL + d_c0 * ue0 * nϵ0 + d_cR * ueR * nϵR
            ∇κ  = d_cL * κL  + d_c0 * κ0  + d_cR * κR
            ∇ϵ  = d_cL * ϵL  + d_c0 * ϵ0  + d_cR * ϵR
            ∇²ϵ = d2_cL * ϵL + d2_c0 * ϵ0 + d2_cR * ϵR

            # Explicit flux term
            F_explicit = flux_factor * ∇nϵue - 10/9 * (κ0 * ∇²ϵ + ∇κ * ∇ϵ)
            
            # Contribution to implicit part from μnϵ * d²ϵ/dz² term
            Aϵ.d[i]    = -10/9 * κ0 * d2_c0 / ne0
            Aϵ.dl[i-1] = -10/9 * κ0 * d2_cL / neL
            Aϵ.du[i]   = -10/9 * κ0 * d2_cR / neR

            # Contribution to implicit part from dμnϵ/dz * ∇ϵ term
            Aϵ.d[i]    -= 10/9 * ∇κ / ne0 * d_c0
            Aϵ.dl[i-1] -= 10/9 * ∇κ / neL * d_cL
            Aϵ.du[i]   -= 10/9 * ∇κ / neR * d_cR

            # Contribution to implicit part from advection term
            Aϵ.d[i]    += flux_factor * ue0 * d_c0
            Aϵ.dl[i-1] += flux_factor * ueL * d_cL
            Aϵ.du[i]   += flux_factor * ueR * d_cR

        if !params.config.LANDMARK && i == 2
            # need to compute heat flux to anode

            mi = params.config.propellant.m
            # 1. current through device
            interior_cell = 3
            je_interior = - e * ne[interior_cell] * params.cache.ue[interior_cell]
            ji_interior = e * sum(Z * U[index.ρiui[Z], interior_cell] for Z in 1:params.config.ncharge) / mi
            jd = ji_interior + je_interior

            # 2. current densities at anode sheath edge
            ji_sheath_edge = e * sum(Z * U[index.ρiui[Z], 1] for Z in 1:params.config.ncharge) / mi
            je_sheath_edge = jd - ji_sheath_edge

            # 3. sheath potential
            Vs = params.cache.ϕ[1] - params.ϕ_L

            # 4. velocity at anode
            ueL = je_sheath_edge / neL / e

            # finite difference derivatives
            #∇nϵue = d_cL * 2/3 * ueL * nϵL + d_c0 * ue0 * nϵ0 + d_cR * ueR * nϵR

            # Explicit flux term
            #F_explicit = flux_factor * ∇nϵue

            # No conduction to anode, compute heat loss due to anode sheath
            Q += d_cL * neL * ueL * Vs

            # Contribution to implicit part from advection term
            # Aϵ.d[i]    = flux_factor * ue0 * d_c0
            # Aϵ.dl[i-1] = flux_factor * 2/3 * ueL * d_cL
            # Aϵ.du[i]   = flux_factor * ueR * d_cR
        end

        # Contribution to implicit part from timestepping
        Aϵ.d[i]    = 1.0 + implicit * Δt * Aϵ.d[i]
        Aϵ.dl[i-1] = implicit * Δt * Aϵ.dl[i-1]
        Aϵ.du[i]   = implicit * Δt * Aϵ.du[i]

        # Explicit right-hand-side
        bϵ[i] = nϵ[i] + Δt * (Q - explicit * F_explicit)
    end

    # Solve equation system using Thomas algorithm
    tridiagonal_solve!(nϵ, Aϵ, bϵ)

    #println("test")
    #@show nϵ[1] / ne[1]
   # @show nϵ[2] / ne[2]

    # Make sure Tev is positive, limit if below user-configured minumum electron temperature
    @inbounds for i in 1:ncells
        if isnan(nϵ[i]) || isinf(nϵ[i]) || nϵ[i] / ne[i] < params.config.min_electron_temperature || nϵ[i] < params.config.min_electron_temperature
            nϵ[i] = 3/2 *  params.config.min_electron_temperature * ne[i]
        end
        Tev[i] = 2/3 * nϵ[i] / ne[i]
        if params.config.LANDMARK
            params.cache.pe[i] = nϵ[i]
        else
            params.cache.pe[i] = 2/3 * nϵ[i]
        end
    end
end