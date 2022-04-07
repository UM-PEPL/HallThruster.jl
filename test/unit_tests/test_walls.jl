@testset "Wall loss tests" begin
    no_losses = HallThruster.NoWallLosses()

    U = nothing
    params = nothing
    i = nothing

    @test no_losses(U, params, i) == 0.0

    landmark_losses = HallThruster.ConstantSheathPotential(-20.0, 1.0, 1.0)

    Tev = 4.0
    ne = 1e18
    cache = (;
        ne = [ne], Tev = [Tev]
    )
    transition_function = HallThruster.StepFunction()
    config = (;
        thruster = (;shielded = true),
        transition_function,
        propellant = HallThruster.Xenon,
    )
    index = (;nϵ = 1)


    params = (;cache, config, index, z_cell = [0.0], L_ch = 0.0)

    U = [ne * Tev * 3/2;;]

    @test no_losses(U, params, 1) == 0.0
    @test landmark_losses(U, params, 1) == 6.0 * 1e7 * exp(-20 / 6.0)

    ce = sqrt(8 * Tev * HallThruster.e / π / HallThruster.me)

    mi = HallThruster.Xenon.m
    me = HallThruster.me
    @test HallThruster.effective_loss_frequency(Tev) == ce / 2

    γ1 = 0.5
    γ2 = 1.0
    @test HallThruster.compute_wall_sheath_potential(Tev, γ1, mi) == -Tev*log(0.5*(1-γ1)*sqrt(2*mi/π/me))
    @test HallThruster.compute_wall_sheath_potential(Tev, γ2, mi) == -1.02 * Tev

    ideal_dielectric = HallThruster.IdealDielectric
    @test HallThruster.SEE_yield(ideal_dielectric, 100.0) == 0.0
    BN = HallThruster.BoronNitride
    @test HallThruster.SEE_yield(BN, 10.0) ≈ BN.Γ * 10.0^BN.b * BN.a

    γ = HallThruster.SEE_yield(BN, Tev)
    Vs = HallThruster.compute_wall_sheath_potential(Tev, γ, mi)
    νe = HallThruster.effective_loss_frequency(Tev)
    sheath_model = HallThruster.WallSheath(BN)

    @test sheath_model(U, params, 1) == νe * Tev * exp(Vs / Tev)

end