Base.@kwdef struct WallMaterial
    name::String
    ϵ_star::Float64
    σ₀::Float64
end

@inline @fastmath function SEE_yield(material::WallMaterial, Tev, γ_max)
    (;σ₀, ϵ_star) = material
    γ = min(γ_max, σ₀ + 1.5 * Tev / ϵ_star * (1 - σ₀))
    return γ
end

const Alumina = WallMaterial(name = "Alumina", ϵ_star = 22, σ₀ = 0.57)
const BoronNitride = WallMaterial(name = "Boron Nitride", ϵ_star = 45, σ₀ = 0.24)
const SiliconDioxide = WallMaterial(name = "Silicon Dioxide", ϵ_star = 18, σ₀ = 0.5)
const BNSiO2 = WallMaterial(name = "Boron Nitride Silicon Dioxide", ϵ_star = 40, σ₀ = 0.54)
const SiliconCarbide = WallMaterial(name = "Silicon Carbide", ϵ_star = 43, σ₀ = 0.69)

function wall_electron_temperature(U, params, i)
    (;cache, config, z_cell) = params

    shielded = config.thruster.shielded

    Tev = cache.Tev[i]

    Tev_channel = shielded * cache.Tev[1] + !shielded * Tev
    Tev_plume = Tev

    L_ch = config.thruster.geometry.channel_length

    Tev = config.transition_function(z_cell[i], L_ch, Tev_channel, Tev_plume)

    return Tev
end

"""
    sheath_potential(Tev, γ, mi))
compute wall sheath to be used for radiative losses and loss to wall.
Goebel Katz equ. 7.3-29, 7.3-44. Assumed nₑuₑ/nᵢuᵢ ≈ 0.5
Sheath potentials are positive by convention in HallThruster.jl.
"""
@inline @fastmath sheath_potential(Tev, γ, mi) = Tev*log((1 - γ) * sqrt(mi/π/me/2))

Base.@kwdef struct WallSheath <: WallLossModel
    material::WallMaterial
    α::Float64 = 0.15
    function WallSheath(material::WallMaterial, α::Float64 = 0.15)
        return new(material, α)
    end
end

function freq_electron_wall(model::WallSheath, U, params, i)
    geometry = params.config.thruster.geometry
    Δr = geometry.outer_radius - geometry.inner_radius
    Tev = wall_electron_temperature(U, params, i)
    mi = params.config.propellant.m

    γ = SEE_yield(model.material, Tev, params.γ_SEE_max)
    params.cache.γ_SEE[i] = γ

    νew = model.α * √(e * Tev / mi) * γ / (Δr * (1 - γ))

    return νew * params.config.transition_function(params.z_cell[i], params.L_ch, 1.0, params.config.electron_plume_loss_scale)
end

function wall_power_loss(model::WallSheath, U, params, i)
    (;config) = params
    mi = config.propellant.m

    Tev = wall_electron_temperature(U, params, i)

    # space charge limited SEE coefficient
    γ = params.cache.γ_SEE[i]

    # Space charge-limited sheath potential
    ϕ_s = sheath_potential(Tev, γ, mi)

    # Compute electron wall collision frequency
    νew = params.cache.νew[i]

    # Compute wall power loss rate
    W = νew * (2 * (1 - 0.5 * model.material.σ₀) * Tev +  (1 - γ) * ϕ_s) / γ

    return W
end
