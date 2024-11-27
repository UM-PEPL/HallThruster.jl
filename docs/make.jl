using Documenter, HallThruster

push!(LOAD_PATH, "../src/")

makedocs(
    sitename = "HallThruster.jl",
    pages = [
        "Home" => "index.md",
        "Release notes" => "NEWS.md",
        "Tutorial" => "run.md",
        "Background" => "background.md",
        "Physics model" => "physics.md",
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
        "Postprocessing" => "postprocessing.md",
        "JSON interface" => "json.md",
        "Grid generation" => "grid.md",
        "Numerics" => "numerics.md",
        "Verification" => "verification.md",
        "Internals" => "internals.md",
        "Contribution" => "contribution.md",
        "Citation" => "citation.md",
    ],
)

deploydocs(
    repo = "github.com/UM-PEPL/HallThruster.jl.git",
)
