
function update_heavy_species!(dU, U, params, t)
    ####################################################################
    #extract some useful stuff from params

    (;index, z_edge, scheme, fluids, species_range_dict, reactions) = params
    (;
        source_neutrals, source_ion_continuity, source_ion_momentum,
        electron_pressure_coupled, min_electron_temperature, propellant
    ) = params.config

    (;ue, μ, F, UL, UR) = params.cache

    ncells = size(U, 2)

    mi = propellant.m

    ##############################################################
    #FLUID MODULE

    #fluid BCs now in update_values struct

    ncharge = params.config.ncharge
    coupled = electron_pressure_coupled
    nϵ = @views U[index.nϵ, :]

    compute_fluxes!(F, UL, UR, U, params)

    @inbounds for i in 2:ncells-1
        left = left_edge(i)
        right = right_edge(i)

        Δz = z_edge[right] - z_edge[left]

        dU[index.ρn, i] = (F[index.ρn, left] - F[index.ρn, right]) / Δz

        # User-provided neutral source term
        dU[index.ρn, i] += source_neutrals(U, params, i)

        ne = 0.0

        for Z in 1:ncharge
            ne += Z * U[index.ρi[Z], i] / mi

            # Acceleration source term
            Q_accel = -Z * e * U[index.ρi[Z], i] / mi * ue[i] / μ[i]

            dU[index.ρi[Z]  , i] = (F[index.ρi[Z],   left] - F[index.ρi[Z],   right]) / Δz
            dU[index.ρiui[Z], i] = (F[index.ρiui[Z], left] - F[index.ρiui[Z], right]) / Δz + Q_accel

            # User-provided source terms
            dU[index.ρi[Z],   i] += source_ion_continuity[Z](U, params, i)
            dU[index.ρiui[Z], i] += source_ion_momentum[Z  ](U, params, i)
        end

        ϵ = max(params.config.min_electron_temperature, U[index.nϵ, i] / ne)

        # Source terms due to ionization
        for r in reactions
            reactant_index = species_range_dict[r.reactant.symbol][1]
            product_index  = species_range_dict[r.product.symbol ][1]
            ρ_reactant = U[reactant_index, i]
            k = r.rate_coeff
            ρdot = k(ϵ) * ρ_reactant * ne
            dU[reactant_index, i] -= ρdot
            dU[product_index, i] += ρdot
        end
    end

    return nothing
end