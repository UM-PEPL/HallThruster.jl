"""
    Thermal_Conductivity
The abstract supertype of all types of thermal conductivity models.
Subtype this to define your own model.
Models return the coefficient in the heat flux `q = -κ * ∇Tev`.
"""
abstract type ThermalConductivityModel end

#=============================================================================
 Serialization
==============================================================================#
@inline function thermal_conductivity_models()
    return (; Braginskii, Mitchner, LANDMARK_conductivity)
end

function Serialization.SType(::Type{T}) where {T <: ThermalConductivityModel}
    return Serialization.TaggedUnion()
end

function Serialization.options(::Type{T}) where {T <: ThermalConductivityModel}
    return thermal_conductivity_models()
end

#=============================================================================
 Built-in model definition
==============================================================================#

"""
    Braginskii conductivity model for fully-ionized plasmasa
    S. I. Braginskii, in Reviews of Plasma Physics (1965), Vol. 1, p. 205.
"""
struct Braginskii <: ThermalConductivityModel end

function (model::Braginskii)(κ, params)
    (; μ, ne, Tev, Z_eff) = params.cache
    @inbounds for i in eachindex(κ)
        # Hunana et al. (2022), Eq. (60): https://arxiv.org/abs/2201.11561
        κ_coef = 13 / 4 + sqrt(2) / Z_eff[i]
        κ[i] = κ_coef * μ[i] * ne[i] * Tev[i]
    end
    return κ
end

"""
    Mitchner-Kruger conductivity model for partially-ionized plasmas
    M. Mitchner and C. H. Kruger, Jr., Partially Ionized Gases (1973). Pg. 94
"""
struct Mitchner <: ThermalConductivityModel end

function (model::Mitchner)(κ, params)
    (; νe, νei, μ, ne, Tev) = params.cache
    @inbounds for i in eachindex(κ)
        # final calculation
        κ[i] = (2.4 / (1 + νei[i] / (√(2) * νe[i]))) * μ[i] * ne[i] * Tev[i]
    end
    return κ
end

"""
    LANDMARK, uses q = -10/9 μnϵ ∇ϵ, where ϵ = 3/2 Tev.
    Expressed as a temperature-gradient coefficient this is 5/2 μneTev.
"""
struct LANDMARK_conductivity <: ThermalConductivityModel end

function (model::LANDMARK_conductivity)(κ, params)
    (; μ, ne, Tev) = params.cache
    @inbounds for i in eachindex(κ)
        κ[i] = (5 / 2) * μ[i] * ne[i] * Tev[i]
    end
    return κ
end
