"""
    Thermal_Conductivity
The abstract supertype of all types of thermal conductivity models.
Subtype this to define your own model.
"""


abstract type ThermalConductivityModel end

"""
Lookup for thermal conductivity coefficients, from Table 1 of 
S. I. Braginskii, in Reviews of Plasma Physics, edited by M. A. Leontovich (Consultants Bureau, New York, 1965), Vol. 1, p. 205.
"""
const LOOKUP_ZS = [1.0, 2.0, 3.0, 4.0, Inf]
const LOOKUP_CONDUCTIVITY_COEFFS = [4.66, 4.0, 3.7, 3.6, 3.2]
const ELECTRON_CONDUCTIVITY_LOOKUP = LinearInterpolation(LOOKUP_ZS, LOOKUP_CONDUCTIVITY_COEFFS)

#=============================================================================
 Begin definition of built-in models
==============================================================================#

"""
    Braginskii, uses closure from   S. I. Braginskii, inReviews of Plasma Physics, edited byM. A. Leontovich (Consultants Bureau, New York, 1965), Vol. 1,p. 205.
    But with a linear interpolation for the charge state
"""
struct Braginskii <: ThermalConductivityModel end

function (model::Braginskii)(κ, params)
    for i in eachindex(κ)
        #get coefficient from charge states
        κ_coef = ELECTRON_CONDUCTIVITY_LOOKUP(params.cache.Z_eff[i])
        #use both classical and anomalous collision frequencies
        ν = (params.cache.νc[i] +  params.cache.νan[i])
        #1/(ωce^2 * me) * q, 1/ term is from Braginskii 1965, q is to convert Te from eV to J
        cyclotron_term = me / (e * (params.cache.B[i])^2)
        #final calculation
        κ[i] = κ_coef * params.cache.ne[i] * params.cache.Tev[i] * ν * cyclotron_term
    end
    return κ
end

"""
    Mitchner, uses closure from M. Mitchner and C. H. Kruger, Jr., Partially Ionized Gases (John Wiley andSons, Inc., New York, 1973). Pg. 94
    Apply factor of 1/(1+Hallparam^2) to convert from parallel to perpendicular direction
"""

struct Mitchner <: ThermalConductivityModel end

function (model::Mitchner)(κ, params)
    for i in eachindex(κ)
        #use both classical and anomalous collision frequencies
        ν = (params.cache.νc[i] +  params.cache.νan[i])
        #calculate mobility using above collision frequency 
        mobility = electron_mobility(ν, params.cache.B[i])
        #final calculation
        κ[i] = (2.4 / (1 + params.cache.νei[i]   / √(2) / ν))   * mobility * params.cache.ne[i] * params.cache.Tev[i]
    end
    return κ
end


"""
    LANDMARK, uses 10/9*
"""

struct LANDMARK_conductivity <: ThermalConductivityModel end

function (model::LANDMARK_conductivity)(κ, params)
    for i in eachindex(κ)
        κ[i] = (10/9) * params.cache.μ[i] * (3/2) * params.cache.ne[i] * params.cache.Tev[i]
    end
    return κ
end