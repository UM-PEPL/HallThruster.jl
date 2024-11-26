using HallThruster: HallThruster as het

function test_json()
    @testset "JSON" begin
        test_path = joinpath(het.PACKAGE_ROOT, "test", "unit_tests")
        json_path = joinpath(test_path, "input_shifted.json")
        sol = het.run_simulation(json_path)
        config = sol.params.config

        @test config.anom_model isa het.LogisticPressureShift
        @test config.anom_model.model.width ≈ 0.0025
        @test config.anom_model.model.center ≈ 0.025
        @test config.anom_model.model.hall_min ≈ 1 / 160
        @test config.anom_model.model.hall_max ≈ 1 / 16
        @test config.anom_model.alpha ≈ 43.0
        @test config.anom_model.pstar ≈ 3e-5
        @test config.anom_model.z0 ≈ -0.003
        @test config.anom_model.dz ≈ 0.005
        @test config.discharge_voltage ≈ 300.0
        @test config.thruster.name == "SPT-100"
        @test config.propellant == het.Xenon
        @test config.anode_mass_flow_rate ≈ 3e-6
        @test config.ion_wall_losses == true
        @test sol.params.simulation.adaptive == true
        @test config.neutral_ingestion_multiplier == 6.0
        @test config.apply_thrust_divergence_correction == false

        json_path = joinpath(test_path, "input_twozone.json")
        sol = het.run_simulation(json_path)
        config = sol.params.config
        @test config.anom_model isa het.TwoZoneBohm
        @test config.anom_model.c1 ≈ 1 / 160
        @test config.anom_model.c2 ≈ 1 / 16
    end
end
