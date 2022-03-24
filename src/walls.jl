abstract type WallLossModel end

struct NoWallLosses <: WallLossModel end

@inline (model::NoWallLosses)(U, params ,i) = 0.0

Base.@kwdef struct ConstantSheathPotential <: WallLossModel
    sheath_potential::Float64 = 20.0
    inner_loss_coeff::Float64
    outer_loss_coeff::Float64
end

function (model::ConstantSheathPotential)(U, params, i)
    (;z_cell, config, index) = params
    L_ch = config.geometry.channel_length
    z = z_cell[i]

    ne = params.cache.ne[i]
    ϵ = U[index.nϵ, i] / ne[i]

    (;sheath_potential, inner_loss_coeff, outer_loss_coeff) = model
    αϵ = params.config.transition_function(z, L_ch, inner_loss_coeff, outer_loss_coeff)
    W = 1e7 * αϵ * ϵ * exp(-sheath_potential / ϵ)

    return W
end

struct WallSheath <: WallLossModel end

@inline function (model::WallSheath)(U, params, i)
    error("Sheath wall loss model not yet implemented.")
end

Base.@kwdef struct WallMaterial
    a::Float64
    b::Float64
    Γ::Float64
end

@inline function SEE_yield(material::WallMaterial, Tev)
    (;a, b, Γ) = material
    return Γ * a * Tev^b
end

const IdealDielectric = WallMaterial(a = 0.0,   b = 0.0,   Γ = 0.0)
const Alumina         = WallMaterial(a = 0.145, b = 0.650, Γ = 1.49)
const BoronNitride    = WallMaterial(a = 0.150, b = 0.549, Γ = 1.38)
const BNSiO2          = WallMaterial(a = 0.123, b = 0.528, Γ = 1.36)
const StainlessSteel  = WallMaterial(a = 0.040, b = 0.610, Γ = 1.44)
