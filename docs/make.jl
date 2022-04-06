using Documenter, HallThruster

push!(LOAD_PATH,"../src/")

makedocs(
    sitename = "HallThruster.jl",
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Physics model" => "physics.md",
        "Configuration" => [
            "config.md",
            "initialization.md",
            "thrusters.md",
            "fluxes.md",
            "collisions.md",
            "anomalous_transport.md",
            "wall_loss_models.md",
            "boundary_conditions.md",
            "source_terms.md",
        ],
        "Numerics" => "numerics.md",
        "Verification" => "verification.md",
        "Internals" => "internals.md"
    ],
)

deploydocs(
    repo = "github.com/UM-PEPL/HallThruster.jl.git",
    devbranch = "main"
)
