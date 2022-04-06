using Documenter, HallThruster

push!(LOAD_PATH,"../src/")

makedocs(
    sitename = "HallThruster.jl",
    pages=[
        "Home" => "index.md",
        "Physics model" => "physics.md",
        "Configuration" => [
            "config.md",
            "initialization.md",
            "thrusters.md",
            "collisions.md",
            "anomalous_transport.md",
            "wall_loss_models.md",
            "boundary_conditions.md",
            "source_terms.md",
        ],
        "Tutorial" => "tutorial.md",
        "Fluxes and Numerics" => "fluxes_and_numerics.md",
        "Verification" => "verification.md",
        "Internals" => "internals.md"
    ],
)

deploydocs(
    repo = "github.com/UM-PEPL/HallThruster.jl.git",
    devbranch = "main"
)
