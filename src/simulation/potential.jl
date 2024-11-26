function compute_electric_field!(∇ϕ, cache, un, apply_drag)
    (; ji, Id, ne, μ, ∇pe, channel_area, ui, νei, νen, νan) = cache

    if (apply_drag)
        (; νei, νen, νan, ui) = cache
    end

    @inbounds for i in eachindex(∇ϕ)
        E = ((Id[] / channel_area[i] - ji[i]) / e / μ[i] - ∇pe[i]) / ne[i]

        if (apply_drag)
            ion_drag = ui[1, i] * (νei[i] + νan[i]) * me / e
            neutral_drag = un * νen[i] * me / e
            E += ion_drag + neutral_drag
        end

        ∇ϕ[i] = -E
    end

    return ∇ϕ
end

function anode_sheath_potential(params)
    (; LANDMARK, anode_boundary_condition) = params.config
    if LANDMARK
        return 0.0
    end

    return anode_sheath_potential(params.cache, anode_boundary_condition)
end

function anode_sheath_potential(cache, anode_bc)
    (; ne, ji, channel_area, Tev, Id) = cache

    # Compute anode sheath potential
    @inbounds if anode_bc == :sheath
        ce = sqrt(8 * e * Tev[1] / π / me)
        je_sheath = e * ne[1] * ce / 4

        # discharge current density
        jd = Id[] / channel_area[1]

        # current densities at sheath edge
        ji_sheath_edge = ji[1]
        je_sheath_edge = jd - ji_sheath_edge

        current_ratio = je_sheath_edge / je_sheath
        if current_ratio ≤ 0.0
            Vs = 0.0
        else
            Vs = -Tev[1] * log(min(1.0, je_sheath_edge / je_sheath))
        end
    else
        Vs = 0.0
    end

    return Vs
end
