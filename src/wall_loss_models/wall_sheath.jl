function wall_electron_temperature(params, i)
    (; cache, config, z_cell) = params

    shielded = config.thruster.shielded

    Tev = cache.Tev[i]

    Tev_channel = shielded * cache.Tev[1] + !shielded * Tev
    Tev_plume = Tev

    L_ch = config.thruster.geometry.channel_length

    Tev = linear_transition(
        z_cell[i], L_ch, config.transition_length, Tev_channel, Tev_plume,)

    return Tev
end

"""
    sheath_potential(Tev, γ, mi))
compute wall sheath to be used for radiative losses and loss to wall.
Goebel Katz equ. 7.3-29, 7.3-44. Assumed nₑuₑ/nᵢuᵢ ≈ 0.5
Sheath potentials are positive by convention in HallThruster.jl.
"""
@inline @fastmath sheath_potential(Tev, γ, mi) = Tev * log((1 - γ) * sqrt(mi / π / me / 2))

Base.@kwdef struct WallSheath <: WallLossModel
    material::WallMaterial
    loss_coefficient::Float64 = 1.0
end

# compute the edge-to-center density ratio
# (c.f https://iopscience.iop.org/article/10.1088/0963-0252/24/2/025017)
# this is approximate, but deviations only become noticable when the knudsen number becomes small, which is not true in our case
@inline edge_to_center_density_ratio() = 0.86 / sqrt(3);

function freq_electron_wall(model::WallSheath, params, i)
    (; config, cache) = params
    (; ncharge) = config
    mi = config.propellant.m
    #compute radii difference
    geometry = config.thruster.geometry
    Δr = geometry.outer_radius - geometry.inner_radius
    #compute electron wall temperature
    Tev = wall_electron_temperature(params, i)
    #calculate and store SEE coefficient
    γ = SEE_yield(model.material, Tev, params.γ_SEE_max)
    cache.see_yield[i] = γ

    # compute the ion current to the walls
    h = edge_to_center_density_ratio()
    j_iw = 0.0
    for Z in 1:ncharge
        niw = h * cache.ni[Z, i]
        j_iw += Z * model.loss_coefficient * niw * sqrt(Z * e * Tev / mi)
    end

    # compute electron wall collision frequency
    νew = j_iw / (Δr * (1 - γ)) / cache.ne[i]

    return νew
end

function wall_power_loss!(Q, ::WallSheath, params)
    (; config, cache, z_cell, L_ch, ncells) = params
    mi = config.propellant.m

    @inbounds for i in 2:(ncells - 1)
        Tev = wall_electron_temperature(params, i)

        # space charge limited SEE coefficient
        γ = params.cache.see_yield[i]

        # Space charge-limited sheath potential
        V_s = sheath_potential(Tev, γ, mi)

        # Compute electron wall collision frequency with transition function for energy wall collisions in plume
        νew = cache.radial_loss_frequency[i] * linear_transition(
            z_cell[i], L_ch, config.transition_length, 1.0, config.electron_plume_loss_scale,
        )

        # Compute wall power loss rate
        Q[i] = νew * (2 * Tev + (1 - γ) * V_s)
    end

    return nothing
end
