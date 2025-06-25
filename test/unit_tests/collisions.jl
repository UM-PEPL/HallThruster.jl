using HallThruster: HallThruster as het

struct Model2 <: het.AnomalousTransportModel end

het.num_anom_variables(::Model2) = 2

function test_collisions()
    Tev = 30 #[eV]
    m = het.Xenon.m #
    Te = Tev * het.e / het.kB #
    ne = 1.0e18 #[#/m^3]
    nn = 0.5e18 #[#/m^3]
    B = 1.0
    ν_an = 0.0
    σ_en = 6.6e-19 * ((Tev / 4 - 0.1) / (1 + (Tev / 4)^1.6)) #[m^2]
    ln_λ = 24 - 0.5 * log(1.0e-6 * ne / Tev^2)
    @test ln_λ ≈ het.coulomb_logarithm(ne, Tev)
    Tev = 9
    ln_λ = 23 - 0.5 * log(1.0e-6 * ne / Tev^3)
    @test ln_λ ≈ het.coulomb_logarithm(ne, Tev)
    ν_c = σ_en * nn * sqrt(8 * het.kB * Te / pi / m) +
        2.9e-12 * ne * ln_λ / (Tev)^1.5
    μ_e = het.e / (het.me * ν_c) /
        (1 + (het.e * B / (het.me * ν_c))^2)
    @test μ_e ≈ het.electron_mobility(ν_an + ν_c, B)

    index = (ρn = [1], ρi = [2], nϵ = 3)
    cache = (;
        nn = [nn], ne = [ne], B = [B], Tev = [Tev], Z_eff = [1.0], νan = [0.0], κ = [0.0],
        μ = μ_e, νc = ν_c,
    )
    c1 = 1 / 160
    c2 = 1 / 16
    anom_model = het.TwoZoneBohm(c1, c2)
    thruster = het.SPT_100
    transition_length = 0.0

    common_opts = (;
        electron_ion_collisions = false, propellant = het.Xenon, anode_mass_flow_rate = 5.0e-6,
        discharge_voltage = 300.0, ncharge = 1, domain = (0.0, 1.0), thruster, transition_length,
    )

    config_landmark = het.Config(;
        anom_model,
        electron_neutral_model = :Landmark,
        common_opts...,
    )
    config_none = het.Config(;
        anom_model,
        electron_neutral_model = :None,
        common_opts...,
    )

    Xe_0 = het.Xenon(0)

    en_landmark = het.load_elastic_collisions(
        config_landmark.electron_neutral_model, [Xe_0],
    )
    en_none = het.load_elastic_collisions(config_none.electron_neutral_model, [Xe_0])

    grid1 = (; cell_centers = [0.02])
    grid2 = (; cell_centers = [0.03])

    params_landmark = (;
        het.params_from_config(config_landmark)...,
        iteration = [0], cache, index, grid = grid1,
        L_ch = thruster.geometry.channel_length, electron_neutral_collisions = en_landmark,
    )
    params_none = (;
        het.params_from_config(config_landmark)...,
        iteration = [0], cache, index, grid = grid2,
        L_ch = thruster.geometry.channel_length, electron_neutral_collisions = en_none,
    )

    @test het.freq_electron_neutral(
        params_landmark.electron_neutral_collisions, nn, Tev,
    ) ≈ 2.5e-13 * nn
    @test het.freq_electron_neutral(params_none.electron_neutral_collisions, nn, Tev) ==
        0.0

    Z = 1
    @test het.freq_electron_ion(ne, Tev, Z) ==
        2.9e-12 * Z^2 * ne * het.coulomb_logarithm(ne, Tev, Z) / Tev^1.5

    config_landmark.anom_model(
        params_landmark.cache.νan, params_landmark, config_landmark,
    )
    config_none.anom_model(params_none.cache.νan, params_none, config_none)

    model = het.NoAnom()

    model(params_landmark.cache.νan, params_landmark, config_landmark)

    @test params_landmark.cache.νan[1] == 0.0

    @test het.ELECTRON_CONDUCTIVITY_LOOKUP(1) == 4.66
    @test het.ELECTRON_CONDUCTIVITY_LOOKUP(1.5) == 4.33

    conductivity_model = het.LANDMARK_conductivity()
    conductivity_model(params_landmark.cache.κ, params_landmark)
    @test params_landmark.cache.κ[1] ≈ 5 / 3 * μ_e * ne * Tev

    @test het.num_anom_variables(model) == 0

    model2 = Model2()
    return @test het.num_anom_variables(model2) == 2
end

test_collisions()
