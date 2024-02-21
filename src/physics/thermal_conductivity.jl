"""
    Thermal_Conductivity
The abstract supertype of all types of thermal conductivity models.
Subtype this to define your own model.
"""


abstract type ThermalConductivityModel end
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
        #use both classical and anomalous collision frequencies
        ν = (params.cache.νc[i] +  params.cache.νan[i])
        #1/(ωce^2 * me) * q, 1/ term is from Braginskii 1965, q is to convert Te from eV to J
        cyclotron_term = me / (e * (params.cache.B[i])^2)

        κ[i] = 4.66 * params.cache.ne[i] * params.cache.Tev[i] * ν * cyclotron_term

    end
    return κ
end