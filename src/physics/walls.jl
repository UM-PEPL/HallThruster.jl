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

    Tev_in = if config.thruster.shielded
        params.cache.Tev[1]
    else
        params.cache.Tev[i]
    end

    Tev_out = params.cache.Tev[i]

    L_ch = config.thruster.geometry.channel_length

    Tev = config.transition_function(z_cell[i], L_ch, Tev_in, Tev_out)

    γ = SEE_yield(model.material, Tev)
    ϕ_s = sheath_potential(Tev, γ, params.config.propellant.m)
    νₑ = effective_loss_frequency(Tev)
    W = νₑ*Tev*exp(ϕ_s/Tev)
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
Space charge limited above γ = 0.99.
"""
function sheath_potential(Tev, γ, mi)
    if γ < 1 - 8.3*sqrt(me/mi)
        # Sheath is not space charge limited
        ϕ_w = -Tev*log(0.5*(1-γ)*sqrt(2*mi/π/me))
    else
        # Sheath is space charge limited
        ϕ_w = -1.02*Tev
    end
    return ϕ_w
end
