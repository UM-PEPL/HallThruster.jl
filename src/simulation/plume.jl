function update_plume_geometry!(U, params; initialize=false)
    (; mi, z_cell, L_ch, exit_plane_index, config, A_ch, cache, index, ncells) = params
    (;
        channel_area, inner_radius, outer_radius, channel_height, dA_dz, tan_div_angle, Tev
    ) = cache

    if !initialize && !params.config.solve_plume
        return nothing
    end

    geometry = config.thruster.geometry
    r_in = geometry.inner_radius
    r_out = geometry.outer_radius

    Tev_exit = Tev[exit_plane_index]

    tan_div_angle[1] = 0.0
    dA_dz[1] = 0.0
    channel_area[1] = A_ch
    inner_radius[1] = r_in
    outer_radius[1] = r_out
    channel_height[1] = r_out - r_in

    for i in 2:ncells
        tanδ_is_zero = z_cell[i] <= L_ch || initialize

        ui =
            sum(U[index.ρiui[Z], i] for Z in 1:(config.ncharge)) /
            sum(U[index.ρi[Z], i] for Z in 1:(config.ncharge))
        tanδ_i = sqrt(5 * e * Tev_exit / 3 / mi) / ui

        tan_div_angle[i] = tanδ_is_zero * 0.0 + !tanδ_is_zero * max(0.0, min(π / 4, tanδ_i))
        avg_tan_δ = 0.5 * (tan_div_angle[i] + tan_div_angle[i - 1])
        Δz = z_cell[i] - z_cell[i - 1]
        inner_radius[i] = max(0.0, inner_radius[i - 1] - avg_tan_δ * Δz)
        outer_radius[i] = outer_radius[i - 1] + avg_tan_δ * Δz
        channel_height[i] = outer_radius[i] - inner_radius[i]
        channel_area[i] = π * (outer_radius[i]^2 - inner_radius[i]^2)
        dA_dz[i] = (channel_area[i] - channel_area[i - 1]) / Δz
    end

    return nothing
end
