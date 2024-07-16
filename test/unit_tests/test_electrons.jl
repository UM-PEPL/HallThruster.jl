@testset "Collisions and mobility" begin
    Tev = 30 #[eV]
    m = HallThruster.Xenon.m #
    Te = Tev*HallThruster.e/HallThruster.kB #
    ne = 1e18 #[#/m^3]
    nn = 0.5e18 #[#/m^3]
    B = 1.0
    ν_an = 0.0
    σ_en = 6.6e-19*((Tev/4 - 0.1)/(1 + (Tev/4)^1.6)) #[m^2]
    @test σ_en ≈ HallThruster.σ_en(Tev)
    ln_λ = 24 - 0.5*log(1e-6*ne/Tev^2)
    @test ln_λ ≈ HallThruster.coulomb_logarithm(ne, Tev)
    Tev = 9
    ln_λ = 23 - 0.5*log(1e-6*ne/Tev^3)
    @test ln_λ ≈ HallThruster.coulomb_logarithm(ne, Tev)
    ν_c = σ_en*nn*sqrt(8*HallThruster.kB*Te/pi/m) + 2.9e-12*ne*ln_λ/(Tev)^1.5
    μ_e = HallThruster.e/(HallThruster.me * ν_c)/(1+(HallThruster.e*B/(HallThruster.me*ν_c))^2)
    @test μ_e ≈ HallThruster.electron_mobility(ν_an + ν_c, B)

    (;e, me) = HallThruster

    mi = HallThruster.Xenon.m

    index = (ρn = [1], ρi = [2], nϵ = 3)
    cache = (;nn = [nn], ne = [ne], B = [B], Tev = [Tev], Z_eff = [1.0], νan = [0.0], κ = [0.0],
             μ = μ_e, νc = ν_c,
    )
    c1 = 1/160
    c2 = 1/16
    anom_model = HallThruster.TwoZoneBohm(c1, c2)
    thruster = HallThruster.SPT_100
    transition_length = 0.0

    config_landmark = (;
        anom_model, propellant = HallThruster.Xenon, electron_neutral_model = HallThruster.LandmarkElectronNeutral(),
        electron_ion_collisions = false, ncharge = 1, thruster, transition_length
    )
    config_gk = (;
        anom_model, propellant = HallThruster.Xenon, electron_neutral_model = HallThruster.GKElectronNeutral(),
        electron_ion_collisions = true, ncharge = 1, thruster, transition_length,
    )
    config_complex = (;
        anom_model, propellant = HallThruster.Xenon, electron_neutral_model = HallThruster.ElectronNeutralLookup(),
        electron_ion_collisions = true, ncharge = 1, thruster, transition_length
    )
    config_none = (;
        anom_model, propellant = HallThruster.Xenon, electron_neutral_model = HallThruster.NoElectronNeutral(),
        electron_ion_collisions = false, ncharge = 1, thruster, transition_length
    )

    Xe_0 = HallThruster.Xenon(0)

    en_landmark = HallThruster._load_reactions(config_landmark.electron_neutral_model, [Xe_0])
    en_gk = HallThruster._load_reactions(config_gk.electron_neutral_model, [Xe_0])
    en_complex = HallThruster._load_reactions(config_complex.electron_neutral_model, [Xe_0])
    en_none = HallThruster._load_reactions(config_none.electron_neutral_model, [Xe_0])

    params_landmark = (;cache, index, config = config_landmark, z_cell = [0.02], L_ch = thruster.geometry.channel_length, electron_neutral_collisions = en_landmark)
    params_gk = (;cache, index, config = config_gk, z_cell = [0.02], L_ch = thruster.geometry.channel_length, electron_neutral_collisions = en_gk)
    params_complex = (;cache, index, config = config_complex, z_cell = [0.03], L_ch = thruster.geometry.channel_length, electron_neutral_collisions = en_complex)
    params_none = (;cache, index, config = config_none, z_cell = [0.03], L_ch = thruster.geometry.channel_length, electron_neutral_collisions = en_none)

    U = [mi * nn; mi * ne; ne * 3/2 * Tev ;;]
    @test HallThruster.freq_electron_neutral(params_landmark, 1) == 2.5e-13 * nn
    @test isapprox(HallThruster.freq_electron_neutral(params_gk, 1), HallThruster.σ_en(Tev) * nn * sqrt(8 * e * Tev / π / me), rtol = 0.01)
    @test HallThruster.freq_electron_neutral(params_none, 1) == 0.0

    Z = 1
    @test HallThruster.freq_electron_ion(ne, Tev, Z) == 2.9e-12 * Z^2 * ne * HallThruster.coulomb_logarithm(ne, Tev, Z) / Tev^1.5
    @test HallThruster.freq_electron_electron(ne, Tev) == 5e-12 * ne * HallThruster.coulomb_logarithm(ne, Tev) / Tev^1.5


    params_landmark.config.anom_model(params_landmark.cache.νan, params_landmark)
    params_none.config.anom_model(params_none.cache.νan, params_none)

    model = HallThruster.NoAnom()

    model(params_landmark.cache.νan, params_landmark)

    @test params_landmark.cache.νan[1] == 0.0

    @test HallThruster.ELECTRON_CONDUCTIVITY_LOOKUP(1) == 4.66
    @test HallThruster.ELECTRON_CONDUCTIVITY_LOOKUP(1.5) == 4.33

    conductivity_model = HallThruster.LANDMARK_conductivity()
    conductivity_model(params_landmark.cache.κ, params_landmark)
    @test params_landmark.cache.κ[1] ≈ 5/3 * μ_e * ne * Tev

    @test HallThruster.num_anom_variables(model) == 0
    @test HallThruster.allocate_anom_variables(model, 2) == Vector{Float64}[]

    struct Model2 <: HallThruster.AnomalousTransportModel end

    HallThruster.num_anom_variables(::Model2) = 2

    model2 = Model2()

    @test HallThruster.num_anom_variables(model2) == 2
    @test HallThruster.allocate_anom_variables(model2, 2) == [
        [0.0, 0.0], [0.0, 0.0]
    ]

    @test HallThruster.allocate_anom_variables(model2, 0) == [
        Float64[], Float64[]
    ]
end
