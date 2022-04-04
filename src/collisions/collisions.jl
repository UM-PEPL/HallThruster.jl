abstract type CollisionModel end

struct NoCollisions <: CollisionModel end

struct SimpleElectronNeutral <: CollisionModel end

struct FullCollisionModel <: CollisionModel end

"""
    freq_electron_neutral(model::CollisionModel, nn, Tev)
Effective frequency of electron scattering caused by collisions with neutrals
"""
@inline freq_electron_neutral(::NoCollisions, nn, Tev)          = 0.0
@inline freq_electron_neutral(::SimpleElectronNeutral, nn, Tev) = 2.5e-13 * nn
@inline freq_electron_neutral(::FullCollisionModel, nn, Tev)    = σ_en(Tev) * nn * sqrt(8 * e * Tev / π / me)

function freq_electron_neutral(U, params, i)
    (;index) = params
    nn = U[index.ρn, i] / params.config.propellant.m
    ne = params.cache.ne[i]
    Tev = params.cache.Tev[i]
    return freq_electron_neutral(params.config.collision_model, nn, Tev)
end

"""
    freq_electron_ion(ne, Tev, Z)
Effective frequency at which electrons are scattered due to collisions with ions
"""
@inline freq_electron_ion(::CollisionModel, ne, Tev, Z)     = 0.0
@inline freq_electron_ion(::FullCollisionModel, ne, Tev, Z) = 2.9e-12 * Z^2 * ne * coulomb_logarithm(ne, Tev, Z) / Tev^1.5

function freq_electron_ion(U, params, i)
    (;index) = params
    mi = params.config.propellant.m
    # Compute effective charge state
    ne = params.cache.ne[i]
    ni_sum = 0.0
    ni_sum = sum(U[index.ρi[Z], i]/mi for Z in 1:params.config.ncharge)
    Z_eff = ne / ni_sum
    Tev = params.cache.Tev[i]
    νei = Tev ≤ 0.0 || ne ≤ 0.0 ? 0.0 : freq_electron_ion(params.config.collision_model, ne, Tev, Z_eff)
    return νei
end

"""
    freq_electron_electron(ne, Tev)
Effective frequency at which electrons are scattered due to collisions with other electrons
"""
@inline freq_electron_electron(ne::Number, Tev::Number) = 5e-12 * ne * coulomb_logarithm(ne, Tev) / Tev^1.5

function freq_electron_electron(U, params, i)
    ne = params.cache.ne[i]
    Tev = params.cache.Tev[i]
    return freq_electron_electron(ne, Tev)
end

freq_electron_wall(U, params, i) = params.config.transition_function(params.z_cell[i], params.L_ch, params.config.wall_collision_freq, 0.0)

freq_electron_anom(U, params, i) = params.config.anom_model(U, params, i)

"""
    coulomb_logarithm(ne, Tev, Z = 1)

calculate coulomb logarithm for electron-ion collisions as a function of ion
charge state Z, electron number density in m^-3, and electron temperature in eV.
"""
@inline function coulomb_logarithm(ne, Tev, Z = 1)
    if Tev < 10 * Z^2
        ln_Λ = 23 - 0.5 * log(1e-6 * ne * Z / Tev^3)
    else
        ln_Λ = 24 - 0.5 * log(1e-6 * ne / Tev^2)
    end

    return ln_Λ
end
