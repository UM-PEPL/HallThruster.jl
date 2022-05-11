@testset "Plume losses" begin
    config = (;thruster = HallThruster.SPT_100, propellant = HallThruster.Krypton, ncharge = 2)
    geom = config.thruster.geometry
    Δr = geom.outer_radius - geom.inner_radius


    Tn = 300.0
    Tev = 3.0

    fluids = [
        HallThruster.Fluid(Xenon(0), HallThruster.ContinuityOnly(u = 300.0, T = Tn)),
        HallThruster.Fluid(Xenon(1), HallThruster.IsothermalEuler(T = Tn)),
        HallThruster.Fluid(Xenon(2), HallThruster.IsothermalEuler(T = Tn)),
    ]

    cache = (;Tev = [Tev, Tev])
    L_ch = geom.channel_length

    z_cell = [L_ch / 2, 2 * L_ch]

    index = (ρn = 1, ρi = [2, 4], ρiui = [3, 5])
    params = (;cache, config, L_ch, fluids, z_cell, index)

    mi = config.propellant.m

    dU = zeros(5, 2)
    u_bohm = sqrt(HallThruster.e * Tev[1] / mi)
    u_thermal_n = sqrt(HallThruster.kB * Tn / 2 / π / mi)

    ρn = 3.0

    ρi_1 = 5.0
    ρi_2 = 1.0
    ui_1 = 10.0
    ui_2 = sqrt(2) * 10.0

    U = [ρn, ρi_1, ρi_1 * ui_1, ρi_2, ρi_2 * ui_2]
    U = hcat(U, U)

    u_bohm_1 = u_bohm
    u_bohm_2 = sqrt(2) * u_bohm

    HallThruster.apply_plume_losses!(dU, U, params, 1)
    HallThruster.apply_plume_losses!(dU, U, params, 2)

    @test dU[1, 1] == 0.0
    @test dU[2, 1] == 0.0
    @test dU[3, 1] == 0.0
    @test dU[4, 1] == 0.0
    @test dU[5, 1] == 0.0

    @test dU[1, 2] == -ρn * u_thermal_n / Δr
    @test dU[2, 2] == -ρi_1 * u_bohm_1 / Δr
    @test dU[3, 2] == -ρi_1 * ui_1 * u_bohm_1 / Δr
    @test dU[4, 2] == -ρi_2 * u_bohm_2 / Δr
    @test dU[5, 2] == -ρi_2 * ui_2 * u_bohm_2 / Δr
end