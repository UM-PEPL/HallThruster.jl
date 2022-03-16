
function update_heavy_species!(dU, U, params, t)
    ####################################################################
    #extract some useful stuff from params

    (;index, z_edge, propellant, scheme, fluids, species_range_dict, reactions) = params
    (;ue, μ, ∇ϕ, ∇pe) = params.cache

    ncells = size(U, 2) - 2

    mi = propellant.m

    ##############################################################
    #FLUID MODULE

    #fluid BCs now in update_values struct

    ncharge = params.config.ncharge

    un, Tn, γn, Rn = fluids[1].conservation_laws.u, fluids[1].conservation_laws.T, γ(fluids[1]), R(fluids[1])
    Ti = fluids[2].conservation_laws.T

    @views ρn = U[index.ρn, :]
    @views nϵ = U[index.nϵ, :]

    coupled = params.config.electron_pressure_coupled

    first_ind = 2
    last_ind = ncells+1

    # Compute heavy species source terms
    @inbounds for i in first_ind:last_ind

        #Compute dU/dt
        left = left_edge(i)
        right = right_edge(i)

        Δz = z_edge[right] - z_edge[left]

        ρn_L = SA[U[index.ρn, i-1]]
        ρn_R = SA[U[index.ρn, i+1]]
        ρn_0 = SA[U[index.ρn, i]]

        # Compute neutral fluxes
        Fn_L = scheme.flux_function(ρn_L, ρn_0, fluids[1], coupled)
        Fn_R = scheme.flux_function(ρn_0, ρn_R, fluids[1], coupled)

        dU[index.ρn, i] = -(Fn_R[1] - Fn_L[1])/Δz

        # Compute electron density
        neL = 0.0
        neR = 0.0
        ne0 = 0.0
        for Z in 1:ncharge
            neL += Z * U[index.ρi[Z], i-1] / mi
            ne0 += Z * U[index.ρi[Z], i] / mi
            neR += Z * U[index.ρi[Z], i+1] / mi
        end

        # Compute electron energy
        ϵL = max(params.config.min_electron_temperature, nϵ[i-1] / neL)
        ϵ0 = max(params.config.min_electron_temperature, nϵ[i] / ne0)
        ϵR = max(params.config.min_electron_temperature, nϵ[i+1] / neR)

        # Compute ion fluxes and source terms
        for Z in 1:ncharge
            UL = SA[U[index.ρi[Z], i-1], U[index.ρiui[Z], i-1]]
            U0 = SA[U[index.ρi[Z], i], U[index.ρiui[Z], i]]
            UR = SA[U[index.ρi[Z], i+1], U[index.ρiui[Z], i+1]]

            fluid = fluids[Z+1]

            FL = scheme.flux_function(UL, U0, fluid, coupled, ϵL, ϵ0, neL, ne0)
            FR = scheme.flux_function(U0, UR, fluid, coupled, ϵ0, ϵR, ne0, neR)

            F_mass, F_momentum = FR[1] - FL[1], FR[2] - FL[2]

            Q_accel = -Z * e * U[index.ρi[Z], i] / mi * ue[i] / μ[i]

            dU[index.ρi[Z], i] = -F_mass / Δz
            dU[index.ρiui[Z], i] = -F_momentum / Δz + Q_accel
        end

        # Source terms due to ionization
        for r in reactions
            reactant_index = species_range_dict[r.reactant.symbol][1]
            product_index = species_range_dict[r.product.symbol][1]
            ρ_reactant = U[reactant_index, i]
            k = r.rate_coeff
            ρdot = k(ϵ0) * ρ_reactant * ne0
            dU[reactant_index, i] -= ρdot
            dU[product_index, i] += ρdot
        end
    end
end