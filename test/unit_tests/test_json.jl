@testset "JSON" begin
    test_path = joinpath(HallThruster.PACKAGE_ROOT, "test", "unit_tests")
    json_path = joinpath(test_path, "input_shifted.json")
    sol = HallThruster.run_simulation(json_path)
    config = sol.params.config

    @test config.anom_model isa HallThruster.ShiftedTwoZoneBohm
    @test config.anom_model.coeffs[1] ≈ 1/160
    @test config.anom_model.coeffs[2] ≈ 1/16
    @test config.anom_model.alpha ≈ 43.0
    @test config.anom_model.pstar ≈ 3e-5
    @test config.anom_model.z0 ≈ -0.003
    @test config.anom_model.dz ≈ 0.005
    @test config.discharge_voltage ≈ 280.0
    @test config.thruster.name == "SPT-100"
    @test config.propellant == Xenon
    @test config.anode_mass_flow_rate ≈ 3e-6
    @test config.solve_background_neutrals == true
    @test config.ion_wall_losses == true
    @test sol.params.adaptive == true

    json_path = joinpath(test_path, "input_twozone.json")
    sol = HallThruster.run_simulation(json_path)
    config = sol.params.config
    @test config.anom_model isa HallThruster.TwoZoneBohm
    @test config.anom_model.coeffs[1] ≈ 1/160
    @test config.anom_model.coeffs[2] ≈ 1/16
end