function wall_electron_temperature(params, transition_length, i)
    (; cache, grid, thruster) = params

    shielded = thruster.shielded

    Tev = cache.Tev[i]

    Tev_channel = shielded * cache.Tev[1] + !shielded * Tev
    Tev_plume = Tev

    L_ch = thruster.geometry.channel_length

    Tev = linear_transition(
        grid.cell_centers[i], L_ch, transition_length, Tev_channel, Tev_plume,
    )

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
    loss_scale::Float64 = 1.0
    function WallSheath(material::WallMaterial, α::Float64 = 0.15)
        return new(material, α)
    end
end

# compute the edge-to-center density ratio
# (c.f https://iopscience.iop.org/article/10.1088/0963-0252/24/2/025017)
# this is approximate, but deviations only become noticable when the knudsen number becomes small, which is not true in our case
@inline edge_to_center_density_ratio() = 0.86 / sqrt(3);

function freq_electron_wall(model::WallSheath, params, i)
    (; cache, ncharge, thruster, mi, transition_length) = params
    #compute radii difference
    geometry = thruster.geometry
    Δr = geometry.outer_radius - geometry.inner_radius
    #compute electron wall temperature
    Tev = wall_electron_temperature(params, transition_length, i)
    #calculate and store SEE coefficient
    γ = SEE_yield(model.material, Tev, params.γ_SEE_max)
    cache.γ_SEE[i] = γ

    # compute the ion current to the walls
    h = edge_to_center_density_ratio()
    j_iw = 0.0
    for Z in 1:ncharge
        niw = h * cache.ni[Z, i]
        j_iw += Z * model.loss_scale * niw * sqrt(Z * e * Tev / mi)
    end

    # compute electron wall collision frequency
    νew = j_iw / (Δr * (1 - γ)) / cache.ne[i]

    return νew
end

function wall_power_loss!(Q, ::WallSheath, params)
    (; cache, grid, thruster, mi, transition_length, plume_loss_scale) = params
    L_ch = thruster.geometry.channel_length

    @inbounds for i in 2:(length(grid.cell_centers) - 1)
        Tev = wall_electron_temperature(params, transition_length, i)

        # space charge limited SEE coefficient
        γ = params.cache.γ_SEE[i]

        # Space charge-limited sheath potential
        ϕ_s = sheath_potential(Tev, γ, mi)

        # Compute electron wall collision frequency with transition function for energy wall collisions in plume
        νew = cache.radial_loss_frequency[i] * linear_transition(
            grid.cell_centers[i], L_ch, transition_length, 1.0, plume_loss_scale,
        )

        # Compute wall power loss rate
        Q[i] = νew * (2 * Tev + (1 - γ) * ϕ_s)
    end

    return nothing
end

function wall_loss_scale(m::WallSheath)
    return m.loss_scale
end
