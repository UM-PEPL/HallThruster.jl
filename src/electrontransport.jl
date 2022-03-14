abstract type AnomalousTransportModel end
abstract type ZeroEquationModel <: AnomalousTransportModel end

struct TwoZoneBohm <: ZeroEquationModel
    coeffs::NTuple{2, Float64}
    TwoZoneBohm(c1, c2) = new((c1, c2))
end

@inline function (model::TwoZoneBohm)(U, params, icell)
    L_ch = params.L_ch
    B = params.cache.B[icell]
    z = params.z_cell[icell]
    c = model.coeffs
    smoothing_length = params.config.smoothing_length

    ωce = e * B / mₑ

    νan = if z < L_ch
        c[1] * ωce + params.νw[1] # +1e7 anomalous wall from Landmark inside channel
    else
        α = smooth_transition(z, L_ch, smoothing_length, params.νw[1]/1e7, params.νw[2]/1e7)
        β = smooth_transition(z, L_ch, smoothing_length, c[1], c[2])
        β * ωce + α * 1e7
    end
    return νan
end

"""
    σ_en(Tev::Float64)

calculation electron neutral collision cross section in m² 
as a function of electron temperature in eV. Eq. 3.6-13, from Fundamentals of 
Electric Propulsion, Goebel and Katz, 2008.

"""
@inline function σ_en(Tev) #Te in eV, from intro to EP, 3.6-13
    return 6.6e-19 * ((Tev / 4 - 0.1) / (1 + (Tev / 4)^1.6)) #[m^2] 
end

"""
    ln_λ(ne::Float64, Tev::Float64)

calculate coulomb logarithm as a function of electron number density and
electron temperature in eV. Eq. 3.6-15, from Fundamentals of 
Electric Propulsion, Goebel and Katz, 2008.

"""
@inline function ln_λ(ne, Tev) #from intro to EP, 3.6-15, or just assume constant 15-25
    return 23 - 0.5 * log(ne * 1e-6 / Tev^3)
end

"""
    get_v_c(Tev::Float64, nn::Float64, ne::Float64, m::Float64)

calculate classical collision frequency, consisting of electron neutral and electron ion collision
frequencies. Eq. 3.6-12 and 3.6-14, from Fundamentals of 
Electric Propulsion, Goebel and Katz, 2008.

"""
@inline function electron_collision_freq(Tev, nn, ne, m) #classical momentum transfer collision frequency, v_ei + v_en, nn neutral number density, ne electron number density, m ion mass
    #v_en = σ_en(Tev) * nn * sqrt(8 * e * Tev / pi / m) # intro to EP, 3.6-12
    #v_ei = 2.9e-12 * ne * ln_λ(ne, Tev) / Tev^1.5 #intro to EP, 3.6-14
    #return v_en + v_ei #2.5e-13*nn #from Hara paper, is similar to formula for v_ei from intro to EP
    return 2.5e-13*nn #from Hara paper and standard in Landmark, is similar to formula for v_ei from intro to EP
end

"""
    B_field(B_max::Float64, z::Float64, L_ch::Float64)

defines magnetic field as a function of position. 
"""
function B_field(B_max, z, L_ch) #same in Landmark and in FFM model Hara
    B = if z < L_ch
        B_max * exp(-0.5 * ((z - L_ch) / (0.011))^2) #for SPT_100
    else
        B_max * exp(-0.5 * ((z - L_ch) / (0.018))^2)
    end
    return B
end

function Te_func(z, L_ch) #will be gone soon
    return 30 * exp(-(2 * (z - L_ch) / 0.033)^2)
end

"""
    cf_electron_transport(v_an::Float64, v_c::Float64, B::Float64)

calculates electron transport according to the generalized Ohm's law
as a function of the classical and anomalous collision frequencies
and the magnetic field.
"""
@inline function electron_mobility(v_an, v_c, B)
    vₑ = v_an + v_c
    Ω = e * B / (mₑ * vₑ)
    return e / (mₑ * vₑ * (1 + Ω^2))
end

"""
    electron_pressure(ne::Float64, Tev::Float64)

ideal gas law for electrons. Tev is in eV.
"""
function electron_pressure(ne, Tev)
    return e * ne * Tev
end