using Documenter, HallThruster

push!(LOAD_PATH, "../src/")

# makedocs(
#     sitename="HallThruster.jl",
#     pages=[
#         "Home" => "index.md",
#         "Release notes" => "NEWS.md",
#         "Tutorial" => "run.md",
#         "Background" => "background.md",
#         "Configuration" => [
#             "config.md",
#             "propellants.md",
#             "initialization.md",
#             "thrusters.md",
#             "fluxes.md",
#             "collisions.md",
#             "anomalous_transport.md",
#             "electron_thermal_conductivity.md",
#             "wall_loss_models.md",
#             "boundary_conditions.md",
#             "source_terms.md",],
#         "Grid generation" => "grid.md",
#         "Postprocessing" => "postprocessing.md",
#         "JSON interface" => "json.md",
#         "Numerics" => "numerics.md",
#         "Physics model" => "physics.md",
#         "Verification" => "verification.md",
#         "Internals" => "internals.md",
#     ],
# )

function section(dir, pairs)
	return [k => joinpath(dir, v) for (k,v) in pairs]
end

makedocs(
	sitename = "HallThruster.jl",
	pages = [
		"Overview" => "index.md",
		"Release notes" => "NEWS.md",
		"Tutorials" => section("tutorials", [
			"Run a simulation" => "simulation.md",
		]),
		"How-to guides" => section("howto", [
			"Use the JSON interface" => "json.md",
			"Run a simulation from python" => "python.md",
			"Add a new propellant" => "new_propellant.md",
			"Implement an anomalous transport model" => "new_anom_model.md",
		]),
		"Explanations" => section("explanation",[
			"Physics model" => "physics.md",
			"Numerics" => "numerics.md",
			"Grid generation" => "grid_generation.md",
			"Timestepping" => "timestepping.md",
			"Validation and verification" => "verification.md",
			"Initialization" => "initialization.md",
			"Quasi-1D plume model" => "plume.md",
		]),
		"Reference" => section("reference", [
			"Configuration" => "config.md",
			"Simulation parameters" => "simulation_options.md",
			"Thrusters" => "thrusters.md",
			"Propellants" => "propellants.md",
			"Anomalous transport models" => "anomalous_transport.md",
			"Wall loss models" => "wall_loss_models.md",
			"Thermal conductivity models" => "electron_thermal_conductivity.md",
			"Collisions and reactions" => "collisions.md",
			"Hyperbolic schemes" => "schemes.md",
			"Postprocessing" => "postprocessing.md",
			"JSON I/O" => "json_interface.md",
			"Internals" => "internals.md",
		]),
	],
	#pagesonly = true,
	format = Documenter.HTML(
		collapselevel = 1,
	)
)

deploydocs(
    repo="github.com/UM-PEPL/HallThruster.jl.git",
)
