abstract type AnomalousTransportModel end
abstract type ZeroEquationModel <: AnomalousTransportModel end

struct NoAnom <: ZeroEquationModel end

@inline (::NoAnom)(U, params, i) = zero(eltype(U))

struct TwoZoneBohm <: ZeroEquationModel
    coeffs::NTuple{2, Float64}
    TwoZoneBohm(c1, c2) = new((c1, c2))
end

@inline function (model::TwoZoneBohm)(U, params, i)
    B = params.cache.B[i]
    ωce = e * B / me

    β = params.config.transition_function(params.z_cell[i], params.L_ch, model.coeffs[1], model.coeffs[2])

    νan = β * ωce
    return νan
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