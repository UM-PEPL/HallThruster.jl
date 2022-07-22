function solve_potential_cell!(ϕ, params)
    (;z_cell, cache, ϕ_L) = params
    (;∇ϕ, Vs) = cache

    cumtrapz!(ϕ, z_cell, ∇ϕ, ϕ_L + Vs[])

    return ϕ
end

function anode_sheath_potential(U, params)
    (;config, index) = params
    (;ne) = params.cache

    mi = params.config.propellant.m

    # Compute anode sheath potential
    if config.LANDMARK
        Vs = 0.0
    else
        ce = sqrt(8 * e * params.cache.Tev[1] / π / me)
        je_sheath = e * ne[1] * ce / 4

        # discharge current density
        jd = params.cache.Id[] / params.A_ch

        # current densities at sheath edge
        ji_sheath_edge = e * sum(Z * U[index.ρiui[Z], 1] for Z in 1:params.config.ncharge) / mi
        je_sheath_edge = jd - ji_sheath_edge

        current_ratio = je_sheath_edge / je_sheath
        if current_ratio ≤ 0.0
            Vs = 0.0
        else
            Vs = -params.cache.Tev[1] * log(min(1.0, je_sheath_edge / je_sheath))
        end
    end

    return Vs
end
