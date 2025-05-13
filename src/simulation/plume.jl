function initialize_plume_geometry(params)
    (; cache, thruster) = params
    (; channel_area, inner_radius, outer_radius, channel_height) = cache
    geometry = thruster.geometry
    r_in = geometry.inner_radius
    r_out = geometry.outer_radius
    A_ch = geometry.channel_area

    @inbounds begin
        @. channel_area = A_ch
        @. inner_radius = r_in
        @. outer_radius = r_out
        @. channel_height = r_out - r_in
    end
end

function update_plume_geometry!(params)
    (; cache, grid, mi, thruster, ncharge) = params
    (; channel_area, inner_radius, outer_radius, channel_height, dA_dz, tanδ, Tev, niui, ni, dlnA_dz) = cache

    L_ch = thruster.geometry.channel_length
    exit_plane_index = max(findfirst(>=(L_ch), grid.cell_centers), 1)
    Tev_exit = Tev[exit_plane_index]
	inv_mi = 1 / mi

    for i in (exit_plane_index + 1):(grid.num_cells)
        ui = sum(niui[Z, i] for Z in 1:ncharge) /
             sum(ni[Z, i] for Z in 1:ncharge)

		thermal_speed = sqrt((5/3) * e * Tev_exit * inv_mi)
		ui = max(ui, thermal_speed)
		tanδ[i] = clamp(thermal_speed / ui, 0.0, 1.0)

        avg_tan_δ = 0.5 * (tanδ[i] + tanδ[i - 1])
        Δz = grid.cell_centers[i] - grid.cell_centers[i - 1]
        inner_radius[i] = max(0.0, inner_radius[i - 1] - avg_tan_δ * Δz)
        outer_radius[i] = outer_radius[i - 1] + avg_tan_δ * Δz
        channel_height[i] = outer_radius[i] - inner_radius[i]
        channel_area[i] = π * (outer_radius[i]^2 - inner_radius[i]^2)
        dA_dz[i] = (channel_area[i] - channel_area[i - 1]) / Δz
        dlnA_dz[i] = dA_dz[i] / channel_area[i]
    end

    return nothing
end
