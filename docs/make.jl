using Documenter, HallThruster

push!(LOAD_PATH, "../src/")

makedocs(
    sitename="HallThruster.jl",
    pages=[
        "Home" => "index.md",
        "Release notes" => "NEWS.md",
        "Tutorial" => "run.md",
        "Background" => "background.md",
        "Configuration" => [
            "config.md",
            "propellants.md",
            "initialization.md",
            "thrusters.md",
            "fluxes.md",
            "collisions.md",
            "anomalous_transport.md",
            "electron_thermal_conductivity.md",
            "wall_loss_models.md",
            "boundary_conditions.md",
            "source_terms.md",],
        "Grid generation" => "grid.md",
        "Postprocessing" => "postprocessing.md",
        "JSON interface" => "json.md",
        "Numerics" => "numerics.md",
        "Physics model" => "physics.md",
        "Verification" => "verification.md",
        "Internals" => "internals.md",
    ],
)

deploydocs(
    repo="github.com/UM-PEPL/HallThruster.jl.git",
)
