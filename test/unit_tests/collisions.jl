using HallThruster: HallThruster as het

struct Model2 <: het.AnomalousTransportModel end

het.num_anom_variables(::Model2) = 2

function test_collisions()
    @testset "Collisions and mobility" begin
        Tev = 30 #[eV]
        m = het.Xenon.m #
        Te = Tev * het.e / het.kB #
        ne = 1e18 #[#/m^3]
        nn = 0.5e18 #[#/m^3]
        B = 1.0
        ν_an = 0.0
        σ_en = 6.6e-19 * ((Tev / 4 - 0.1) / (1 + (Tev / 4)^1.6)) #[m^2]
        ln_λ = 24 - 0.5 * log(1e-6 * ne / Tev^2)
        @test ln_λ ≈ het.coulomb_logarithm(ne, Tev)
        Tev = 9
        ln_λ = 23 - 0.5 * log(1e-6 * ne / Tev^3)
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

        config_landmark = (;
            anom_model, propellant = het.Xenon,
            electron_neutral_model = :Landmark,
            electron_ion_collisions = false, ncharge = 1, thruster, transition_length,
        )
        config_none = (;
            anom_model, propellant = het.Xenon,
            electron_neutral_model = :None,
            electron_ion_collisions = false, ncharge = 1, thruster, transition_length,
        )

        Xe_0 = het.Xenon(0)

        en_landmark = het.load_elastic_collisions(
            config_landmark.electron_neutral_model, [Xe_0],)
        en_none = het.load_elastic_collisions(config_none.electron_neutral_model, [Xe_0])

        grid1 = (; cell_centers = [0.02])
        grid2 = (; cell_centers = [0.03])

        params_landmark = (;
            iteration = [0], cache, index, config = config_landmark, grid = grid1,
            L_ch = thruster.geometry.channel_length, electron_neutral_collisions = en_landmark,)
        params_none = (;
            iteration = [0], cache, index, config = config_none, grid = grid2,
            L_ch = thruster.geometry.channel_length, electron_neutral_collisions = en_none,)

        @test het.freq_electron_neutral(params_landmark, 1) ≈ 2.5e-13 * nn
        @test het.freq_electron_neutral(params_none, 1) == 0.0

        Z = 1
        @test het.freq_electron_ion(ne, Tev, Z) ==
              2.9e-12 * Z^2 * ne * het.coulomb_logarithm(ne, Tev, Z) / Tev^1.5
        @test het.freq_electron_electron(ne, Tev) ==
              5e-12 * ne * het.coulomb_logarithm(ne, Tev) / Tev^1.5

        params_landmark.config.anom_model(params_landmark.cache.νan, params_landmark)
        params_none.config.anom_model(params_none.cache.νan, params_none)

        model = het.NoAnom()

        model(params_landmark.cache.νan, params_landmark)

        @test params_landmark.cache.νan[1] == 0.0

        @test het.ELECTRON_CONDUCTIVITY_LOOKUP(1) == 4.66
        @test het.ELECTRON_CONDUCTIVITY_LOOKUP(1.5) == 4.33

        conductivity_model = het.LANDMARK_conductivity()
        conductivity_model(params_landmark.cache.κ, params_landmark)
        @test params_landmark.cache.κ[1] ≈ 5 / 3 * μ_e * ne * Tev

        @test het.num_anom_variables(model) == 0
        @test het.allocate_anom_variables(model, 2) == Vector{Float64}[]

        model2 = Model2()

        @test het.num_anom_variables(model2) == 2
        @test het.allocate_anom_variables(model2, 2) == [
            [0.0, 0.0], [0.0, 0.0],
        ]

        @test het.allocate_anom_variables(model2, 0) == [
            Float64[], Float64[],
        ]
    end
end
