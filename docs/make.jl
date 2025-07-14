using HallThruster
using Documenter
using DocumenterCitations

push!(LOAD_PATH, "../src/")

function section(dir, pairs)
    return [k => joinpath(dir, v) for (k, v) in pairs]
end
bib = CitationBibliography(
    joinpath(@__DIR__, "src", "assets", "refs.bib");
    style = :numeric,
)

makedocs(
    sitename = "HallThruster.jl",
    pages = [
        "Overview" => "index.md",
        "Release notes" => "NEWS.md",
        "Science highlights" => "highlights.md",
        "Tutorials" => section(
            "tutorials", [
                "Run a simulation" => "simulation.md",
            ]
        ),
        "How-to guides" => section(
            "howto", [
                "Use JSON input and output" => "json.md",
                "Run a simulation from python" => "python.md",
                "Add a new propellant" => "new_propellant.md",
                "Implement an anomalous transport model" => "new_anom_model.md",
            ]
        ),
        "Explanations" => section(
            "explanation", [
                "Physics model" => "physics.md",
                "Numerics" => "numerics.md",
                "Grid generation" => "grid_generation.md",
                "Timestepping" => "timestepping.md",
                "Validation and verification" => "verification.md",
                "Initialization" => "initialization.md",
                "Quasi-1D plume model" => "plume.md",
            ]
        ),
        "Reference" => section(
            "reference", [
                "Configuration" => "config.md",
                "Simulation parameters" => "simparams.md",
                "Solution" => "solution.md",
                "Postprocessing" => "postprocessing.md",
                "Thrusters" => "thrusters.md",
                "Propellants" => "propellants.md",
                "Anomalous transport models" => "anomalous_transport.md",
                "Wall loss models" => "wall_loss_models.md",
                "Thermal conductivity models" => "electron_thermal_conductivity.md",
                "Collisions and reactions" => "collisions.md",
                "Internals" => "internals.md",
            ]
        ),
    ],
    pagesonly = true,
    format = Documenter.HTML(
        collapselevel = 1,
        assets = [
            "assets/citation.css",
            "assets/favicon.ico",
        ]
    ),
    modules = [HallThruster],
    checkdocs = :exports,
    plugins = [bib],
)

deploydocs(
    repo = "github.com/UM-PEPL/HallThruster.jl.git",
)
