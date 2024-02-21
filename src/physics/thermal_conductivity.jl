"""
    Thermal_Conductivity
The abstract supertype of all types of thermal conductivity models.
Subtype this to define your own model.
"""


abstract type ThermalConductivityModel end

"""
Lookup for thermal conductivity coefficients, from Table 1 of 
S. I. Braginskii, inReviews of Plasma Physics, edited byM. A. Leontovich (Consultants Bureau, New York, 1965), Vol. 1,p. 205.
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