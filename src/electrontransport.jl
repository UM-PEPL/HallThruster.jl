abstract type AnomalousTransportModel end
abstract type ZeroEquationModel <: AnomalousTransportModel end

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

struct DataDriven <: ZeroEquationModel
    coeffs::NTuple{1, Float64}
    DataDriven(c1) = new((c1,))
end

@inline function (model::DataDriven)(U, params, icell)
    (;index) = params
    (;∇ϕ, B, νan) = params.cache
    c = model.coeffs

    ui = abs(U[index.ρiui[1], icell] / U[index.ρi[1], icell])
    ωce = e * B[icell] / me
    vde = max(ui, abs(-∇ϕ[icell] / B[icell]))
    if νan[icell] == 0.0
        α = 1.0
    else
        α = 0.5
    end
    return α * max(1e-4 * ωce, c[1] * ωce * ui / vde) + (1-α) * νan[icell]
end

"""
    electron_mobility(νan::Float64, νc::Float64, B::Float64)

calculates electron transport according to the generalized Ohm's law
as a function of the classical and anomalous collision frequencies
and the magnetic field.
"""
@inline function electron_mobility(νan, νc, B)
    νe = νan + νc
    Ω = e * B / (me * νe)
    return e / (me * νe * (1 + Ω^2))
end
