using HallThruster: HallThruster as het
using JSON3: JSON3

include(joinpath(het.TEST_DIR, "unit_tests", "serialization_test_utils.jl"))

function test_json()
    @testset "JSON" begin
        test_path = joinpath(het.PACKAGE_ROOT, "test", "json")
        json_path = joinpath(test_path, "input_shifted.json")
        sol = het.run_simulation(json_path)
        config = sol.config

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
        config = sol.config
        @test config.anom_model isa het.TwoZoneBohm
        @test config.anom_model.c1 ≈ 1 / 160
        @test config.anom_model.c2 ≈ 1 / 16

        outfile = "output.json"
        @test ispath(outfile)

        out = JSON3.read(outfile)
        @test haskey(out, "input")
        @test haskey(out.input, "config")
        @test haskey(out.input, "simulation")
        @test haskey(out.input, "postprocess")
        @test haskey(out, "output")
        @test haskey(out.output, "error")
        @test haskey(out.output, "retcode")
        @test haskey(out.output, "average")
        @test haskey(out.output, "frames")

        # Test that reading the output file produces the same inputs we originally ran the simulation with
        new_sol = het.run_simulation(outfile)
        @test struct_eq(new_sol.config, sol.config)
        @test struct_eq(new_sol.params.simulation, sol.params.simulation)
        @test struct_eq(new_sol.params.postprocess, sol.params.postprocess)

        # Test restarting simulation from a json file
        # restart should be a different solution than the original,
        # since it runs for an additional `duration`
        restart = het.run_simulation(outfile, restart = outfile)
        @test !isapprox(new_sol.savevals[end].Id[], restart.savevals[end].Id[])

        rm(outfile, force = true)
    end
end
