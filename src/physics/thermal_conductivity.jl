"""
    Thermal_Conductivity
The abstract supertype of all types of thermal conductivity models.
Subtype this to define your own model.
"""
abstract type ThermalConductivityModel end

"""
Lookup for thermal conductivity coefficients, from Table 1 of
S. I. Braginskii, in Reviews of Plasma Physics (1965), Vol. 1, p. 205.
"""
const LOOKUP_ZS = range(1, 5, length = 5)
const LOOKUP_CONDUCTIVITY_COEFFS = [4.66, 4.0, 3.7, 3.6, 3.2]
const ELECTRON_CONDUCTIVITY_LOOKUP = LinearInterpolation(
    LOOKUP_ZS, LOOKUP_CONDUCTIVITY_COEFFS,
)

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
struct Braginskii <: ThermalConductivityModel
    anom_scale::Float64
    function Braginskii(anom_scale::Float64)
        @check_in_interval anom_scale 0 1
        return new(anom_scale)
    end
end

Braginskii(;anom_scale::Float64=1.0) = Braginskii(anom_scale)

function (model::Braginskii)(κ, params)
    (; νan, νc, ne, Tev, Z_eff, B) = params.cache
    @inbounds for i in eachindex(κ)
        # get coefficient from charge states
        κ_coef = ELECTRON_CONDUCTIVITY_LOOKUP(Z_eff[i])
        ωce = e * B[i] / me
        ν_eff = νc[i] + model.anom_scale * νan[i]
        Ω = ωce / ν_eff
        μ = (e / me) / (ν_eff * (1 + Ω^2))
        # final calculation
        κ[i] = κ_coef *  ne[i] * Tev[i] * μ
    end
    return κ
end

"""
    Mitchner-Kruger conductivity model for partially-ionized plasmas
    M. Mitchner and C. H. Kruger, Jr., Partially Ionized Gases (1973). Pg. 94
"""
struct Mitchner <: ThermalConductivityModel
    anom_scale::Float64
    function Mitchner(anom_scale::Float64)
        @check_in_interval anom_scale 0 1
        return new(anom_scale)
    end
end

Mitchner(;anom_scale::Float64=1.0) = Mitchner(anom_scale)


function (model::Mitchner)(κ, params)
    (; νei, νan, νc, ne, Tev, B) = params.cache
    @inbounds for i in eachindex(κ)
        ωce = e * B[i] / me
        ν_eff = νc[i] + model.anom_scale * νan[i]
        Ω = ωce / ν_eff
        μ = e / (me * ν_eff * (1 + Ω^2))
        κ[i] = (2.4 / (1 + νei[i] / (√(2) * ν_eff))) * μ * ne[i] * Tev[i]
    end
    return κ
end

"""
    LANDMARK, uses 10/9 μnϵ
"""
struct LANDMARK_conductivity <: ThermalConductivityModel end

function (model::LANDMARK_conductivity)(κ, params)
    (; μ, ne, Tev) = params.cache
    @inbounds for i in eachindex(κ)
        κ[i] = (10 / 9) * μ[i] * (3 / 2) * ne[i] * Tev[i]
    end
    return κ
end
