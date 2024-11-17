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
    LOOKUP_ZS, LOOKUP_CONDUCTIVITY_COEFFS)

#=============================================================================
 Begin definition of built-in models
==============================================================================#

"""
    Braginskii conductivity model for fully-ionized plasmasa
    S. I. Braginskii, in Reviews of Plasma Physics (1965), Vol. 1, p. 205.
"""
struct Braginskii <: ThermalConductivityModel end

function (model::Braginskii)(κ, params)
    (; nu_class, nu_wall, nu_anom, B, ne, Tev, Z_eff) = params.cache
    @inbounds for i in eachindex(κ)
        #get coefficient from charge states
        κ_coef = ELECTRON_CONDUCTIVITY_LOOKUP(Z_eff[i])
        #use both classical and anomalous collision frequencies
        ν = nu_class[i] + nu_wall[i] + nu_anom[i]
        # calculate mobility using above collision frequency
        μ = electron_mobility(ν, B[i])
        #final calculation
        κ[i] = κ_coef * μ * ne[i] * Tev[i]
    end
    return κ
end

"""
    Mitchner-Kruger conductivity model for partially-ionized plasmas
    M. Mitchner and C. H. Kruger, Jr., Partially Ionized Gases (1973). Pg. 94
"""
struct Mitchner <: ThermalConductivityModel end

function (model::Mitchner)(κ, params)
    (; nu_class, nu_wall, nu_ei, nu_anom, B, ne, Tev) = params.cache
    @inbounds for i in eachindex(κ)
        # use both classical and anomalous collision frequencies
        ν = nu_class[i] + nu_wall[i] + nu_anom[i]
        # calculate mobility using above collision frequency
        μ = electron_mobility(ν, B[i])
        # final calculation
        κ[i] = (2.4 / (1 + nu_ei[i] / (√(2) * ν))) * μ * ne[i] * Tev[i]
    end
    return κ
end

"""
    LANDMARK, uses 10/9 μnϵ
"""
struct LANDMARK_conductivity <: ThermalConductivityModel end

function (model::LANDMARK_conductivity)(κ, params)
    (; mobility, ne, Tev) = params.cache
    @inbounds for i in eachindex(κ)
        κ[i] = (10 / 9) * mobility[i] * (3 / 2) * ne[i] * Tev[i]
    end
    return κ
end
