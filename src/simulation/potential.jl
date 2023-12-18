function solve_potential!(ϕ, params)
    (;z_cell, cache, ϕ_L) = params
    (;∇ϕ, Vs) = cache

    cumtrapz!(ϕ, z_cell, ∇ϕ, ϕ_L + Vs[])

    return ϕ
end

function anode_sheath_potential(U, params)
    (;config, index) = params
    (;ne, channel_area) = params.cache

    mi = config.propellant.m
    Ti = config.ion_temperature

    if params.config.LANDMARK
        return 0.0
    end

    ci = sqrt(2 * kB * Ti / mi)

    # Compute anode sheath potential
    if config.anode_boundary_condition == :sheath
        ce = sqrt(8 * e * params.cache.Tev[1] / π / me)
        je_sheath = e * ne[1] * ce / 4

        # discharge current density
        jd = params.cache.Id[] / channel_area[1]

        # current densities at sheath edge
        ji_sheath_edge = e * sum(Z * U[index.ρiui[Z], 1] for Z in 1:params.config.ncharge) / mi

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
