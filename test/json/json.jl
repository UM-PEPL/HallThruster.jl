using JSON3: JSON3
using HallThruster: HallThruster as het

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
@test haskey(out.output.average, "ni")
@test haskey(out.output.average, "niui")
@test haskey(out.output.average, "ui")
@test haskey(out.output.average, "nn")
@test haskey(out.output.average, "ions")
@test haskey(out.output.average, "neutrals")

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

out = JSON3.read(outfile);
output = out.output
@test haskey(output.average, "ions")
@test haskey(output.average, "neutrals")
@test !haskey(output.average, "nn")
@test !haskey(output.average, "ni")
@test !haskey(output.average, "niui")
@test !haskey(output.average, "ui")
for prop in config.propellants
    prop_str = string(prop.gas.short_name)
    @test haskey(output.average.ions, prop_str)
    @test haskey(output.average.neutrals, prop_str)
    @test length(output.average.ions[prop_str]) == prop.max_charge
end

rm(outfile, force = true)
