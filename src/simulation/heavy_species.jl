
function update_heavy_species!(dU, U, params, t)
    ####################################################################
    #extract some useful stuff from params

    (;index, z_edge, config, cache) = params
    (;
        source_neutrals, source_ion_continuity, source_ion_momentum,
        propellant, scheme, thruster, ncharge, ion_wall_losses
    ) = config

    (;F, UL, UR) = cache

    mi = propellant.m


    ncells = size(U, 2)
    nedges = length(z_edge)

    Δr = thruster.geometry.outer_radius - thruster.geometry.inner_radius

    ##############################################################
    #FLUID MODULE

    #fluid BCs now in update_values struct

    compute_fluxes!(F, UL, UR, U, params)

    if scheme.WENO #with WENO 5 only going up to 2nd last to boundary cells
        @inbounds for i in 1:nedges #only fluxes ncells - 1

            # If stencil goes past boundary, we have ghost edges
            # Where the boundary flux is
            i_minus_2 = max(1, i - 2)
            i_minus_1 = max(1, i - 1)
            i_plus_1 = min(nedges, i+1)
            i_plus_2 = min(nedges, i+2)

            # Handle neutrals
            for j in 1:params.num_neutral_fluids
                F[index.ρn[j], i] = WENO5_compute_fluxes(
                F[index.ρn[j], i_minus_2], F[index.ρn[j], i_minus_1], F[index.ρn[j], i], F[index.ρn[j], i_plus_1], F[index.ρn[j], i_plus_2])
            end

            #Handle ions
            for Z in 1:ncharge
                @views @. F[index.ρi[Z]:index.ρiui[Z], i] = WENO5_compute_fluxes(
                    SA[F[index.ρi[Z], i_minus_2], F[index.ρiui[Z], i_minus_2]],
                    SA[F[index.ρi[Z], i_minus_1], F[index.ρiui[Z], i_minus_1]],
                    SA[F[index.ρi[Z], i],   F[index.ρiui[Z], i]],
                    SA[F[index.ρi[Z], i_plus_1], F[index.ρiui[Z], i_plus_1]],
                    SA[F[index.ρi[Z], i_plus_2], F[index.ρiui[Z], i_plus_2]],
                )
            end
        end
    end

    @inbounds for i in 2:ncells-1
        left = left_edge(i)
        right = right_edge(i)

        Δz = z_edge[right] - z_edge[left]

        # Handle neutrals
        for j in 1:params.num_neutral_fluids
            # Neutral fluxes
            dU[index.ρn[j], i] = (F[index.ρn[j], left] - F[index.ρn[j], right]) / Δz

            # User-provided neutral source term
            dU[index.ρn[j], i] += source_neutrals[j](U, params, i)
        end

        # Handle ions
        for Z in 1:ncharge
            # Ion fluxes
            dU[index.ρi[Z]  , i] = (F[index.ρi[Z],   left] - F[index.ρi[Z],   right]) / Δz
            dU[index.ρiui[Z], i] = (F[index.ρiui[Z], left] - F[index.ρiui[Z], right]) / Δz

            # User-provided ion source terms
            dU[index.ρi[Z],   i] += source_ion_continuity[Z](U, params, i)
            dU[index.ρiui[Z], i] += source_ion_momentum[Z  ](U, params, i)
        end

        apply_ion_acceleration!(dU, U, params, i)
        apply_reactions!(dU, U, params, i)

        if ion_wall_losses
            apply_ion_wall_losses!(dU, U, params, i)
        end

        dU[index.nϵ, i] = 0.0
    end

    @. @views dU[:, 1] = 0.0
    @. @views dU[:, end] = 0.0

    return nothing
end
