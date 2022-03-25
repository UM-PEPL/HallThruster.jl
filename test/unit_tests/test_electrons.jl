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
    @test μ_e ≈ HallThruster.electron_mobility(ν_an, ν_c, B)

    (;e, me) = HallThruster

    mi = HallThruster.Xenon.m

    index = (ρn = 1, ρi = [2], nϵ = 3)
    cache = (;ne = [ne], B = [B])
    anom_model = HallThruster.TwoZoneBohm(1/160, 1/16)
    thruster = HallThruster.SPT_100
    transition_function = HallThruster.StepFunction()
    config_simple = (;anom_model, propellant = HallThruster.Xenon, collision_model = HallThruster.SimpleElectronNeutral(), ncharge = 1, thruster, transition_function)
    config_complex = (;anom_model, propellant = HallThruster.Xenon, collision_model = HallThruster.FullCollisionModel(), ncharge = 1, thruster, transition_function)
    config_none = (;anom_model, propellant = HallThruster.Xenon, collision_model = HallThruster.NoCollisions(), ncharge = 1, thruster, transition_function)

    params_1 = (;cache, index, config = config_simple, z_cell = [0.02], L_ch = thruster.geometry.channel_length)
    params_2 = (;cache, index, config = config_complex, z_cell = [0.03], L_ch = thruster.geometry.channel_length)
    params_3 = (;cache, index, config = config_none, z_cell = [0.03], L_ch = thruster.geometry.channel_length)

    U = [mi * nn; mi * ne; ne * 3/2 * Tev ;;]
    @test HallThruster.freq_electron_neutral(U, params_1, 1) == 2.5e-13 * nn
    @test HallThruster.freq_electron_neutral(U, params_2, 1) == HallThruster.σ_en(Tev) * nn * sqrt(8 * e * Tev / π / me)
    @test HallThruster.freq_electron_neutral(U, params_3, 1) == 0.0

    Z = 1

    @test HallThruster.freq_electron_ion(U, params_1, 1) == 0.0
    @test HallThruster.freq_electron_ion(U, params_2, 1) == 2.9e-12 * Z^2 * ne * HallThruster.coulomb_logarithm(ne, Tev, Z) / Tev^1.5
    @test HallThruster.freq_electron_ion(U, params_3, 1) == 0.0

    @test HallThruster.freq_electron_electron(ne, Tev) == 5e-12 * ne * HallThruster.coulomb_logarithm(ne, Tev) / Tev^1.5
    @test HallThruster.freq_electron_electron(U, params_1, 1) == 5e-12 * ne * HallThruster.coulomb_logarithm(ne, Tev) / Tev^1.5

    @test HallThruster.freq_electron_anom(U, params_1, 1) == e/me * 1/160
    @test HallThruster.freq_electron_anom(U, params_2, 1) == e/me * 1/16

    model = HallThruster.NoAnom()
    @test model(U, params_1, 1) == 0.0
    @test model(U, params_1, 1) == 0.0

end