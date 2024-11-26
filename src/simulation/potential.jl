function compute_electric_field!(∇ϕ, cache, un, apply_drag)
    (; ji, Id, ne, μ, ∇pe, channel_area, ui, νei, νen, νan) = cache

    if (apply_drag)
        (; νei, νen, νan, ui) = cache
    end

    for i in eachindex(∇ϕ)
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
    (; config) = params
    (; ne, niui, channel_area) = params.cache

    if params.config.LANDMARK
        return 0.0
    end

    # Compute anode sheath potential
    if config.anode_boundary_condition == :sheath
        ce = sqrt(8 * e * params.cache.Tev[1] / π / me)
        je_sheath = e * ne[1] * ce / 4

        # discharge current density
        jd = params.cache.Id[] / channel_area[1]

        # current densities at sheath edge
        ji_sheath_edge = e * sum(Z * niui[Z, 1] for Z in 1:(params.config.ncharge))
        je_sheath_edge = jd - ji_sheath_edge

        current_ratio = je_sheath_edge / je_sheath
        if current_ratio ≤ 0.0
            Vs = 0.0
        else
            Vs = -params.cache.Tev[1] * log(min(1.0, je_sheath_edge / je_sheath))
        end
    else
        Vs = 0.0
    end

    return Vs
end
