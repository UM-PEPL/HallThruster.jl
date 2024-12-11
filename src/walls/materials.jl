@kwdef struct WallMaterial
    name::String
    ϵ_star::Float64
    σ₀::Float64
end

@inline @fastmath function SEE_yield(material::WallMaterial, Tev, γ_max)
    (; σ₀, ϵ_star) = material
    γ = min(γ_max, σ₀ + 1.5 * Tev / ϵ_star * (1 - σ₀))
    return γ
end

const Alumina = WallMaterial(name = "Alumina", ϵ_star = 22, σ₀ = 0.57)
const BoronNitride = WallMaterial(name = "Boron Nitride", ϵ_star = 45, σ₀ = 0.24)
const SiliconDioxide = WallMaterial(name = "Silicon Dioxide", ϵ_star = 18, σ₀ = 0.5)
const BNSiO2 = WallMaterial(name = "Boron Nitride Silicon Dioxide", ϵ_star = 40, σ₀ = 0.54)
const SiliconCarbide = WallMaterial(name = "Silicon Carbide", ϵ_star = 43, σ₀ = 0.69)

#=============================================================================
 Serialization
==============================================================================#
const wall_materials = (; Alumina, BoronNitride, SiliconDioxide, BNSiO2, SiliconCarbide)
Serialization.SType(::Type{WallMaterial}) = Serialization.Enum()
Serialization.options(::Type{WallMaterial}) = wall_materials
