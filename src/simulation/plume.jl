function initialize_plume_geometry!(params)
    geometry = params.config.thruster.geometry
    r_in = geometry.inner_radius
    r_out = geometry.outer_radius
    A_ch = geometry.channel_area

    (; channel_area, inner_radius, outer_radius, channel_height) = params.cache

    @inbounds begin
        @. channel_area = A_ch
        @. inner_radius = r_in
        @. outer_radius = r_out
        @. channel_height = r_out - r_in
    end
    return
end

function update_plume_geometry!(params)
    (; z_cell, config, cache, ncells) = params
    (; channel_area, inner_radius, outer_radius, channel_height, dA_dz, tan_div_angle, Tev) = cache

    if config.solve_plume
        return
    end

    mi = config.propellant.m
    geometry = config.thruster.geometry
    L_ch = geometry.channel_length

    # TODO: look this code over and see if we can compute this in the loop below
    exit_plane_index = max(findfirst(>(L_ch), z_cell) - 1, 2)
    Tev_exit = Tev[exit_plane_index]

    # Compute divergence angle and geometry
    for i in exit_plane_index:ncells
        # TODO: update this. the logic is kind of convoluted, but essentially,
        # we should first compute set tan_delta to zero to initialize everything (1st time calling)
        # then, we should only loop over the cells after the exit plane
        # The initialization logic should be removed entirely or made simpler,
        # i.e. inner_radius .= geometry.inner_radius somewhere else, not in this loop

        # get charge-weighted axial ion velocity
        ui = sum(cache.niui[Z, i] for Z in 1:(config.ncharge)) /
             sum(cache.ni[Z, i] for Z in 1:(config.ncharge))

        tanδ_i = sqrt(5 * e * Tev_exit / 3 / mi) / ui
        tan_div_angle[i] = max(0.0, tanδ_i)
        avg_tan_δ = 0.5 * (tan_div_angle[i] + tan_div_angle[i - 1])

        Δz = z_cell[i] - z_cell[i - 1]
        increment = avg_tan_δ * Δz

        inner_radius[i] = max(0.0, inner_radius[i - 1] - increment)
        outer_radius[i] = outer_radius[i - 1] + increment

        channel_height[i] = outer_radius[i] - inner_radius[i]
        channel_area[i] = π * (outer_radius[i]^2 - inner_radius[i]^2)

        dA_dz[i] = (channel_area[i] - channel_area[i - 1]) / Δz
    end

    return nothing
end
