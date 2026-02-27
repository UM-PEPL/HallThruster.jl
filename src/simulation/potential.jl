function update_electric_field!(∇ϕ, cache, apply_drag)
    (; ji, Id, ne, μ, ∇pe, channel_area, νei, νen, νan, avg_ion_vel, avg_neutral_vel) = cache

    @inbounds for i in eachindex(∇ϕ)
        E = ((Id[] / channel_area[i] - ji[i]) / e / μ[i] - ∇pe[i]) / ne[i]

        if (apply_drag)
            for s in ion_species
                ion_drag += u_s[i] * νei_s[i] * me / e
            end
            neutral_drag = avg_neutral_vel[i] * νen[i] * me / e
            E += ion_drag + neutral_drag
        end

        ∇ϕ[i] = -E
    end

    return ∇ϕ
end

function integrate_potential!(ϕ, ∇ϕ, grid, V_L)
    # We need to make sure the potential is integrated from the left edge to the right edge,
    # rather than from the left ghost cell center to the right ghost cell center

    # Temorarily replace left and right electric field and grid values with edge values
    EL, ER = ∇ϕ[1], ∇ϕ[end]
    zL, zR = grid.cell_centers[1], grid.cell_centers[end]

    grid.cell_centers[1] = grid.edges[1]
    grid.cell_centers[end] = grid.edges[end]

    ∇ϕ[1] = 0.5 * (EL + ∇ϕ[2])
    ∇ϕ[end] = 0.5 * (ER + ∇ϕ[end - 1])

    # Integrate potential from left to right edge
    cumtrapz!(ϕ, grid.cell_centers, ∇ϕ, V_L)

    # Extrapolate potential to ghost cells
    ϕ[1] = ϕ[1] + (ϕ[1] - ϕ[2])
    ϕ[end] = ϕ[end] + (ϕ[end] - ϕ[end - 1])

    # Replace electric field and cell center values
    grid.cell_centers[1], grid.cell_centers[end] = zL, zR
    ∇ϕ[1], ∇ϕ[end] = EL, ER
    return
end


function anode_sheath_potential(params)
    if params.landmark
        return 0.0
    end
    (; anode_bc, cache) = params
    (; ne, ji, channel_area, Tev, Id) = cache

    # Compute anode sheath potential
    @inbounds if anode_bc == :sheath

        Te_sheath_edge = 0.5 * (Tev[1] + Tev[2])
        ne_sheath_edge = 0.5 * (ne[1] + ne[2])
        ce = sqrt(8 * e * Te_sheath_edge / π / me)
        je_sheath = e * ne_sheath_edge * ce / 4

        # discharge current density
        jd = Id[] / channel_area[1]

        # current densities at sheath edge
        ji_sheath_edge = 0.5 * (ji[1] + ji[2])
        je_sheath_edge = jd - ji_sheath_edge

        current_ratio = je_sheath_edge / je_sheath
        if current_ratio ≤ 0.0
            Vs = 0.0
        else
            Vs = -Te_sheath_edge * log(min(1.0, je_sheath_edge / je_sheath))
        end
    else
        Vs = 0.0
    end

    return Vs
end
