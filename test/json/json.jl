using HallThruster: HallThruster as het, JSON

include("$(het.TEST_DIR)/unit_tests/serialization_test_utils.jl")

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
@test config.anom_model.pstar ≈ 3.0e-5
@test config.anom_model.z0 ≈ -0.12
@test config.anom_model.dz ≈ 0.2
@test config.discharge_voltage ≈ 300.0
@test config.thruster.name == "SPT-100"
@test config.propellants[1].gas == het.Xenon
@test config.propellants[1].flow_rate_kg_s ≈ 3.0e-6
@test config.ion_wall_losses == true
@test sol.simulation.adaptive == true
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

out = JSON.parsefile(outfile)
@test haskey(out, "input")
input = out["input"]
@test haskey(input, "config")
@test haskey(input, "simulation")
@test haskey(input, "postprocess")
@test haskey(out, "output")
output = out["output"]
@test haskey(output, "error")
@test haskey(output, "retcode")
@test haskey(output, "average")
@test haskey(output, "frames")
avg = output["average"]
@test haskey(avg, "ni")
@test haskey(avg, "niui")
@test haskey(avg, "ui")
@test haskey(avg, "nn")
@test haskey(avg, "ions")
@test haskey(avg, "neutrals")

# Test that reading the output file produces the same inputs we originally ran the simulation with
new_sol = het.run_simulation(outfile)
@test struct_eq(new_sol.config, sol.config)
@test struct_eq(new_sol.simulation, sol.simulation)
@test struct_eq(new_sol.postprocess, sol.postprocess)

# Test restarting simulation from a json file
# restart should be a different solution than the original,
# since it runs for an additional `duration`
restart = het.run_simulation(outfile, restart = outfile)
@test !isapprox(new_sol.frames[end].discharge_current[], restart.frames[end].discharge_current[])

#==============================================================================
    Multiple propellants
==============================================================================#
json_path = joinpath(test_path, "input_multiprop.json")
sol = het.run_simulation(json_path)
config = sol.config
@test length(config.propellants) == 2
prop1 = config.propellants[1]
prop2 = config.propellants[2]
@test prop1.gas == het.Xenon
@test prop2.gas == het.Krypton

outfile = "output.json"
@test ispath(outfile)

out = JSON.parsefile(outfile);
avg = out["output"]["average"]
@test haskey(avg, "ions")
@test haskey(avg, "neutrals")
@test !haskey(avg, "nn")
@test !haskey(avg, "ni")
@test !haskey(avg, "niui")
@test !haskey(avg, "ui")
for prop in config.propellants
    prop_str = string(prop.gas.short_name)
    @test haskey(avg.ions, prop_str)
    @test haskey(avg.neutrals, prop_str)
    @test length(avg.ions[prop_str]) == length(prop.allowed_charges)
end

rm(outfile, force = true)
