"""
$(TYPEDEF)

Abstract type for wall loss models in the electron energy equation. Types included with HallThruster.jl are:

- `NoWallLosses`
    No electron energy losses to the wall.
- `ConstantSheathPotential(sheath_potential, inner_loss_coeff, outer_loss_coeff)`
    Power to the walls scales with `α * exp(-sheath_potential))`, where `α = inner_loss_coeff` inside the channel and `α = outer_loss_coeff` outside.
- `WallSheath(material::WallMaterial)`
    Power to the walls is computed self-consistently using a sheath model, dependent on the secondary electron emission yield of the provided material.
    See `WallMaterial` for material options.
"""
abstract type WallLossModel end

struct NoWallLosses <: WallLossModel end

@inline (model::NoWallLosses)(U, params ,i) = 0.0

Base.@kwdef struct ConstantSheathPotential <: WallLossModel
    sheath_potential::Float64 = -20.0
    inner_loss_coeff::Float64
    outer_loss_coeff::Float64
end

function (model::ConstantSheathPotential)(U, params, i)
    (;z_cell, index, L_ch) = params

    z = z_cell[i]

    ne = params.cache.ne[i]
    ϵ = U[index.nϵ, i] / ne

    (;sheath_potential, inner_loss_coeff, outer_loss_coeff) = model
    αϵ = params.config.transition_function(z, L_ch, inner_loss_coeff, outer_loss_coeff)
    W = 1e7 * αϵ * ϵ * exp(sheath_potential / ϵ)

    return W
end

"""
$(TYPEDEF)

Struct containing secondary electron emission (SEE) yield fit coefficients for a material. Used in the `WallSheath` wall loss model

# Fields

$(TYPEDFIELDS)

The SEE yield is computed as a function of electron temperature in eV (`Tev`) as `Γ * a * Tev^b`.

# Available materials

- `IdealDielectric` (a material with zero SEE)
- `Alumina`
- `BoronNitride`
- `BNSiO2`
- `StainlessSteel`

"""
Base.@kwdef struct WallMaterial
    a::Float64
    b::Float64
    Γ::Float64
end

"""
    SEE_yield(material::WallMaterial, Tev)
fit function for SEE with different wall materials
Goebel Katz equ. 7.3-30
"""
@inline function SEE_yield(material::WallMaterial, Tev)
    (;a, b, Γ) = material
    return Γ * a * Tev^b
end

const IdealDielectric = WallMaterial(a = 0.0,   b = 0.0,   Γ = 0.0)
const Alumina         = WallMaterial(a = 0.145, b = 0.650, Γ = 1.49)
const BoronNitride    = WallMaterial(a = 0.150, b = 0.549, Γ = 1.38)
const BNSiO2          = WallMaterial(a = 0.123, b = 0.528, Γ = 1.36)
const StainlessSteel  = WallMaterial(a = 0.040, b = 0.610, Γ = 1.44)

Base.@kwdef struct WallSheath <: WallLossModel
    material::WallMaterial
end

function (model::WallSheath)(U, params, i)
    (;config, z_cell) = params
    mi = config.propellant.m
    Δr = config.thruster.geometry.outer_radius - config.thruster.geometry.inner_radius

    Tev_channel = if config.thruster.shielded
        params.cache.Tev[1]
    else
        params.cache.Tev[i]
    end

    Tev_plume = params.cache.Tev[i]

    L_ch = config.thruster.geometry.channel_length

    α_channel = 1.0
    α_plume = 0.0

    in_channel = config.transition_function(z_cell[i], L_ch, 1.0, 0.0)
    in_plume = 1 - in_channel
    Tev = in_channel * Tev_channel + in_plume * Tev_plume
    α = in_channel * α_channel + in_plume * α_plume
    γ = SEE_yield(model.material, Tev)
    ϕ_s = sheath_potential(Tev, γ, mi)

    #νₑ = effective_loss_frequency(Tev)
    #W = νₑ*Tev*exp(ϕ_s/Tev)

    u_bohm = sqrt(e * params.cache.Z_eff[i] * Tev / mi)
    νew = 2 * u_bohm / Δr / (1 - γ)
    W = α * νew * (2Tev + (1 - γ) * ϕ_s)

    return W
end

"""
    effective_loss_frequency(Tev)
compute effective loss frequency for electron power loss to walls. 
based on Half Maxwellian approximation and correction factor 2 for 
energy lost per electron compared to average energy. Goebel-Katz
equ. 7.3 - 45. 
"""
function effective_loss_frequency(Tev)
    νₑ = 1/4*sqrt(8*e*Tev/π/me)*2
    return νₑ
end

"""
    sheath_potential(Tev, γ, mi))
compute wall sheath to be used for radiative losses and loss to wall.
Goebel Katz equ. 7.3-29, 7.3-44. Assumed nₑuₑ/nᵢuᵢ ≈ 0.5
Space charge limited when γ > 1 - 8.3 √(mᵢ/mₑ).
"""
function sheath_potential(Tev, γ, mi)

    sqrt_me_mi = sqrt(me/mi)

    # space charge limited SEE coefficient
    γ_limited = min(γ, 1 - 8.3 * sqrt_me_mi)

    # space charge-limited sheath potential.
    # by convention in HallThruster.jl, sheath potentials are positive when ion-attracting
    ϕ_w = Tev*log((1 - γ_limited) * sqrt(mi/π/me/2))

    return ϕ_w
end
