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

"""
    σ_en(Tev)

Electron neutral collision cross section in m² as a function of electron temperature in eV.
Eq. 3.6-13, from Fundamentals of Electric Propulsion, Goebel and Katz, 2008.

"""
@inline function σ_en(Tev)
    return 6.6e-19 * ((Tev / 4 - 0.1) / (1 + (Tev / 4)^1.6))
end

"""
    freq_electron_neutral(nn, Tev)
Effective frequency of electron scattering caused by collisions with neutrals
"""
@inline function freq_electron_neutral(nn::Number, Tev::Number, model)
    if model == :simple
        return 2.5e-13 * nn
    else
        σ_en(Tev) * nn * sqrt(8 * e * Tev / π / mₑ)
    end
end

function freq_electron_neutral(U, params, i)
    (;index) = params
    nn = U[index.ρn, i] / params.config.propellant.m
    ne = params.cache.ne[i]
    Tev = 2/3 * U[index.nϵ, i] / ne
    return freq_electron_neutral(nn, Tev, params.config.collision_model)
end

"""
    freq_electron_ion(ne, Tev, Z)
Effective frequency at which electrons are scattered due to collisions with ions
"""
@inline freq_electron_ion(ne::Number, Tev::Number, Z::Number) = 2.9e-12 * Z^2 * ne * coulomb_logarithm(ne, Tev, Z) / Tev^1.5

function freq_electron_ion(U, params, i)
    (;index) = params
    mi = params.config.propellant.m
    # Compute effective charge state
    ne = params.cache.ne[i]
    ni_sum = 0.0
    ni_sum = sum(U[index.ρi[Z], i]/mi for Z in 1:params.config.ncharge)
    Z_eff = ne / ni_sum

    Tev = 2/3 * U[index.nϵ, i] / ne
    return freq_electron_ion(ne, Tev, Z_eff)
end

"""
    freq_electron_electron(ne, Tev)
Effective frequency at which electrons are scattered due to collisions with other electrons
"""
@inline freq_electron_electron(ne::Number, Tev::Number) = 5e-12 * ne * coulomb_logarithm(ne, Tev) / Tev^1.5

function freq_electron_electron(U, params, i)
    ne = params.cache.ne[i]
    Tev = 2/3 * U[params.index.nϵ, i] / ne
    return freq_electron_electron(ne, Tev)
end

freq_electron_wall(U, params, i) = 1e7 * params.config.transition_function(params.z_cell[i], params.L_ch, params.αw[1], 0.0)

freq_electron_anom(U, params, i) = params.anom_model(U, params, i)