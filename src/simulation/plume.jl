function initialize_plume_geometry(params)
    (; cache, grid, thruster, transition_length) = params
    (; channel_area, inner_radius, outer_radius, channel_height, wall_transition) = cache
    geometry = thruster.geometry
    r_in = geometry.inner_radius
    r_out = geometry.outer_radius
    A_ch = geometry.channel_area
    L_ch = geometry.channel_length

    @inbounds begin
        @. channel_area = A_ch
        @. inner_radius = r_in
        @. outer_radius = r_out
        @. channel_height = r_out - r_in
        fill!(cache.dA_dz, 0.0)
        fill!(cache.dlnA_dz, 0.0)
        fill!(cache.tanδ, 0.0)
        cache.channel_area_old .= channel_area
        cache.plume_energy_old .= cache.nϵ
        for i in eachindex(wall_transition)
            wall_transition[i] = linear_transition(
                grid.cell_centers[i], L_ch, transition_length, 1.0, 0.0,
            )
        end

        # Cache the final interior cell affected by ion wall losses. The wall
        # transition and grid are fixed throughout the simulation.
        last_wall_cell = length(wall_transition) - 1
        while last_wall_cell >= 2 && iszero(wall_transition[last_wall_cell])
            last_wall_cell -= 1
        end
        params.last_wall_cell = last_wall_cell
    end

    return nothing
end

@inline function plume_expansion_slope(Tev, mass, velocity)
    thermal_speed = sqrt((5 / 3) * e * Tev / mass)
    return iszero(thermal_speed) ? 0.0 : thermal_speed / max(velocity, thermal_speed)
end

function update_plume_geometry!(params)
    (; cache, grid, thruster) = params
    (; channel_area, inner_radius, outer_radius, channel_height, dA_dz, tanδ, Tev, m_eff, dlnA_dz) = cache

    L_ch = thruster.geometry.channel_length
    geometry = thruster.geometry
    z = grid.cell_centers
    # Sample the reference state at the physical exit, independently of how
    # the exit happens to align with the cell centers.
    Tev_exit = interpolate(L_ch, z, Tev)
    mass_exit = interpolate(L_ch, z, m_eff)
    velocity_exit = interpolate(L_ch, z, cache.avg_ion_vel)
    previous_slope = plume_expansion_slope(Tev_exit, mass_exit, velocity_exit)
    previous_z = L_ch
    spread = 0.0

    @inbounds for i in eachindex(z)
        if z[i] <= L_ch
            inner_radius[i] = geometry.inner_radius
            outer_radius[i] = geometry.outer_radius
            channel_height[i] = geometry.outer_radius - geometry.inner_radius
            channel_area[i] = geometry.channel_area
            tanδ[i] = 0.0
            dA_dz[i] = 0.0
            dlnA_dz[i] = 0.0
            continue
        end

        # Integrate the existing variable-angle closure from the actual exit.
        slope = plume_expansion_slope(Tev_exit, m_eff[i], cache.avg_ion_vel[i])
        spread += 0.5 * (slope + previous_slope) * (z[i] - previous_z)
        tanδ[i] = slope
        inner_radius[i] = max(0.0, geometry.inner_radius - spread)
        outer_radius[i] = geometry.outer_radius + spread
        channel_height[i] = outer_radius[i] - inner_radius[i]
        channel_area[i] = π * (outer_radius[i]^2 - inner_radius[i]^2)
        # dr_out/dz = slope; dr_in/dz = -slope until the hole closes.
        # Once r_in = 0, its contribution to dA/dz is also zero.
        dA_dz[i] = 2π * slope * (outer_radius[i] + inner_radius[i])
        dlnA_dz[i] = dA_dz[i] / channel_area[i]
        previous_slope = slope
        previous_z = z[i]
    end

    return nothing
end
