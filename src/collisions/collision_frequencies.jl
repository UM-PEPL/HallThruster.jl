"""
    freq_electron_neutral(coll::ElasticCollision, nn, Tev)
Effective frequency of electron scattering caused by collisions with neutrals
"""
@inline function freq_electron_neutral(coll::ElasticCollision, nn::Number, Tev::Number)
    return rate_coeff(coll, 3 / 2 * Tev) * nn
end


"""
    freq_electron_neutral!(νen, coll::ElasticCollision, neutral::Fluid, Tev)
Effective frequency of electron scattering caused by collisions with neutrals of a specific species
"""
function freq_electron_neutral!(νen::Vector{T}, coll::ElasticCollision, fluid::FluidContainer, Tev::Vector{T}) where {T <: Number}
    inv_m = inv(fluid.species.element.m)
    @inbounds for i in interior_cells(νen)
        νen[i] += freq_electron_neutral(coll, fluid.density[i] * inv_m, Tev[i])
    end
    νen[1] = νen[2]
    νen[end] = νen[end - 1]
    return νen
end

"""
    freq_electron_ion(ni, Tev, Z, ne = abs(Z) * ni)
Effective electron scattering frequency for ions of number density `ni` and
charge `Z`. `ne` sets electron screening in the Coulomb logarithm.
"""
@inline function freq_electron_ion(ni::Number, Tev::Number, Z::Number, ne::Number = abs(Z) * ni)
    charge = abs(float(Z))
    return 2.9e-12 * charge^2 * ni * coulomb_logarithm(ne, Tev, charge) / sqrt(Tev^3)
end

function freq_electron_ion!(
        νei::Vector{T}, ne::Vector{T}, Tev::Vector{T}, Z::Vector{T},
    ) where {T <: Number}
    @inbounds for i in interior_cells(νei)
        νei[i] = freq_electron_ion(ne[i] / abs(Z[i]), Tev[i], Z[i], ne[i])
    end
    νei[1] = νei[2]
    νei[end] = νei[end - 1]
    return
end

function freq_electron_ion!(νei, ne, Tev, fluids::AbstractVector{FluidContainer})
    fill!(νei, 0.0)
    @inbounds for fluid in fluids
        Z = fluid.species.Z
        iszero(Z) && continue
        inv_m = inv(fluid.species.element.m)
        for i in interior_cells(νei)
            ni = fluid.density[i] * inv_m
            νei[i] += freq_electron_ion(ni, Tev[i], Z, ne[i])
        end
    end
    νei[1] = νei[2]
    νei[end] = νei[end - 1]
    return nothing
end

"""
    $(TYPEDSIGNATURES)
Update the classical collision frequency.
"""
function freq_electron_classical!(
        νc::Vector{T}, νen::Vector{T}, νei::Vector{T},
        νiz::Vector{T}, νex::Vector{T}, landmark::Bool,
    ) where {T}
    @inbounds for i in eachindex(νc)
        νc[i] = νen[i] + νei[i]
    end

    if landmark
        return
    end

    @inbounds for i in eachindex(νc)
        νc[i] += νiz[i] + νex[i]
    end

    return
end

"""
    coulomb_logarithm(ne, Tev, Z = 1)

calculate coulomb logarithm for electron-ion collisions as a function of ion
charge state Z, electron number density in m^-3, and electron temperature in eV.
"""
@inline function coulomb_logarithm(ne, Tev, Z = 1)
    if Tev < 10 * Z^2
        ln_Λ = 23 - 0.5 * log(1.0e-6 * abs(ne) * Z^2 / Tev^3)
    else
        ln_Λ = 24 - 0.5 * log(1.0e-6 * abs(ne) / Tev^2)
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
