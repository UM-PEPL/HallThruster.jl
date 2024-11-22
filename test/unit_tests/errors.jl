using HallThruster: HallThruster as het

@testset "Errors and Warnings" begin
    Landmark_config = het.Config(;
        thruster = het.SPT_100,
        domain = (0.0u"cm", 8.0u"cm"),
        discharge_voltage = 300.0u"V",
        anode_mass_flow_rate = 5u"mg/s",
        LANDMARK = true,
    )

    @test_throws ErrorException het.run_simulation(
        Landmark_config; dt = 5e-9, duration = 4e-9,
        grid = het.EvenGrid(2), nsave = 10,)

    config = het.Config(;
        thruster = het.SPT_100,
        domain = (0.0u"cm", 8.0u"cm"),
        discharge_voltage = 300.0u"V",
        anode_mass_flow_rate = 5u"mg/s",
    )

    @test_logs (:warn,
        "CFL for adaptive timestepping set higher than stability limit of 0.8. Setting CFL to 0.799.",) het.run_simulation(
        config; dt = 5e-9, duration = 0e-9, grid = het.EvenGrid(2),
        nsave = 10, adaptive = true, CFL = 0.9,)
end
