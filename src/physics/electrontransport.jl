abstract type AnomalousTransportModel end
abstract type ZeroEquationModel <: AnomalousTransportModel end

struct NoAnom <: ZeroEquationModel end

@inline (::NoAnom)(U, params, i) = zero(eltype(U))

struct TwoZoneBohm <: ZeroEquationModel
    coeffs::NTuple{2, Float64}
    TwoZoneBohm(c1, c2) = new((c1, c2))
end

struct JornsModel <: ZeroEquationModel
    coeffs::NTuple{3, Float64}
    JornsModel(c1, c2, c3) = new((c1, c2, c3))
end


@inline function (model::TwoZoneBohm)(U, params, i)
    B = params.cache.B[i]
    ωce = e * B / me

    β = params.config.transition_function(params.z_cell[i], params.L_ch, model.coeffs[1], model.coeffs[2])

    νan = β * ωce
    return νan
end

@inline function (model::JornsModel)(U, params, i)
    index = params.index
    B = params.cache.B[i]
    ωce = e * B / me
    E = -params.cache.∇ϕ[i]
    c_s = params.cache.Tev[i]
    u_i = U[index.ρiui[1], i]/U[index.ρi[1], i]

    νan = ωce*(model.coeffs[1] + model.coeffs[2]*abs(u_i)/(model.coeffs[3]*c_s + abs(E/B)))
    return νan
end


struct DataDriven <: HallThruster.ZeroEquationModel
    coeffs::Vector{Float64}
end

@inline function (model::DataDriven)(U, params, i)
    (;index) = params
    (;∇ϕ, B, νan, Tev, νc, ue, νan) = params.cache
    c = model.coeffs

    ncells = size(U, 2)

    mi = params.config.propellant.m
    ui = abs(U[index.ρiui[1], i] / U[index.ρi[1], i])
    ωce = HallThruster.e * B[i] / HallThruster.me
    cs = sqrt(HallThruster.e * Tev[i] / mi)

    #Ωₑ = ωce / (νc[i] + νan[i])
    #vde = abs(Ωₑ * ue[i])
    vde = abs(∇ϕ[i] / B[i])

    f_an = max(0.0, ωce * (c[2] * ui / (c[3] * cs + vde)))

    α = 0.5
    # smooth anomalous transport in time and space
    if i == 1
        f_an = α * f_an + (1 - α) * (νan[1] + νan[2])/2
    elseif i == ncells
        f_an = α * f_an + (1 - α) * (νan[end] + νan[end-1])/2
    else
        f_an = α * f_an + (1 - α) * (νan[i-1] + 2 * νan[i] + νan[i+1])/4
    end
    return f_an
end



"""
    electron_mobility(νan::Float64, νc::Float64, B::Float64)

calculates electron transport according to the generalized Ohm's law
as a function of the classical and anomalous collision frequencies
and the magnetic field.
"""
@inline function electron_mobility(νan, νc, B)
    νe = νan + νc
    return electron_mobility(νe, B)
end

@inline function electron_mobility(νe, B)
    Ω = e * B / (me * νe)
    return e / (me * νe * (1 + Ω^2))
end


@inline electron_sound_speed(Tev) = sqrt(8 * e * Tev / π / me)