function solve_potential!(ϕ, params)
    (; z_cell, cache, V_L) = params
    (; electric_field, Vs) = cache

    cumtrapz!(ϕ, z_cell, electric_field, V_L + Vs[])

    return ϕ
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
