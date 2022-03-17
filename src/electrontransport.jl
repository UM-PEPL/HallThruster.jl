abstract type AnomalousTransportModel end
abstract type ZeroEquationModel <: AnomalousTransportModel end

struct TwoZoneBohm <: ZeroEquationModel
    coeffs::NTuple{2, Float64}
    TwoZoneBohm(c1, c2) = new((c1, c2))
end

@inline function (model::TwoZoneBohm)(U, params, i)
    B = params.cache.B[i]
    ωce = e * B / mₑ

    β = params.config.transition_function(params.z_cell[i], params.L_ch, model.coeffs[1], model.coeffs[2])

    νan = β * ωce
    return νan
end

struct DataDriven <: ZeroEquationModel
    coeffs::NTuple{1, Float64}
    DataDriven(c1) = new((c1,))
end


@inline function (model::DataDriven)(U, params, icell)
    (;index) = params
    (;∇ϕ, B, νan) = params.cache
    c = model.coeffs

    ui = abs(U[index.ρiui[1], icell] / U[index.ρi[1], icell])
    ωce = e * B[icell] / mₑ
    vde = max(ui, abs(-∇ϕ[icell] / B[icell]))
    if νan[icell] == 0.0
        α = 1.0
    else
        α = 0.5
    end
    return α * max(1e-4 * ωce, c[1] * ωce * ui / vde) + (1-α) * νan[icell]
end


"""
    electron_collision_freq(Tev::Float64, nn::Float64, ne::Float64, m::Float64)

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