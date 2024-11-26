"""
    freq_electron_neutral(model::ElectronNeutralModel, nn, Tev)
Effective frequency of electron scattering caused by collisions with neutrals
"""
@inline function freq_electron_neutral(
        collisions::Vector{ElasticCollision}, nn::Number, Tev::Number,)
    νen = 0.0
    @inbounds for c in collisions
        νen += rate_coeff(c, 3 / 2 * Tev) * nn
    end
    return νen
end

function freq_electron_neutral!(
        νen::Vector{T}, collisions::Vector{ElasticCollision},
        nn::Vector{T}, Tev::Vector{T},) where {T <: Number}
    @inbounds for i in eachindex(νen)
        νen[i] = freq_electron_neutral(collisions, nn[i], Tev[i])
    end
end

"""
    freq_electron_ion(ne, Tev, Z)
Effective frequency at which electrons are scattered due to collisions with ions
"""
@inline function freq_electron_ion(ne::T, Tev::T, Z::T) where {T <: Number}
    return 2.9e-12 * Z^2 * ne * coulomb_logarithm(ne, Tev, Z) / sqrt(Tev^3)
end

function freq_electron_ion!(
        νei::Vector{T}, ne::Vector{T}, Tev::Vector{T}, Z::Vector{T},) where {T <: Number}
    @inbounds for i in eachindex(νei)
        νei[i] = freq_electron_ion(ne[i], Tev[i], Z[i])
    end
end

"""
    $(TYPEDSIGNATURES)
Update the classical collision frequency.
"""
function freq_electron_classical!(
        νc::Vector{T}, νen::Vector{T}, νei::Vector{T},
        νiz::Vector{T}, νex::Vector{T}, landmark::Bool,) where {T}
    @inbounds for i in eachindex(νc)
        νc[i] = νen[i] + νei[i]
    end

    if landmark
        return
    end

    @inbounds for i in eachindex(νc)
        νc[i] += νiz[i] + νex[i]
    end
end

"""
    coulomb_logarithm(ne, Tev, Z = 1)

calculate coulomb logarithm for electron-ion collisions as a function of ion
charge state Z, electron number density in m^-3, and electron temperature in eV.
"""
@inline function coulomb_logarithm(ne, Tev, Z = 1)
    if Tev < 10 * Z^2
        ln_Λ = 23 - 0.5 * log(1e-6 * ne * Z^2 / Tev^3)
    else
        ln_Λ = 24 - 0.5 * log(1e-6 * ne / Tev^2)
    end

    return ln_Λ
end

"""
    electron_mobility(νe, B)

calculates electron transport according to the generalized Ohm's law
as a function of sum of the classical and anomalous collision frequencies
and the magnetic field.
"""
@inline function electron_mobility(νe, B)
    Ω = e * B / (me * νe)
    return e / (me * νe * (1 + Ω^2))
end

@inline electron_sound_speed(Tev) = sqrt(8 * e * Tev / π / me)
