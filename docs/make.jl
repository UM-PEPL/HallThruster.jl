using Documenter, HallThruster

push!(LOAD_PATH,"../src/")

makedocs(
    sitename = "HallThruster.jl",
    pages=[
        "Home" => "index.md",
        "Configuration" => Any[
            "config.md",
            "initialization.md",
            "collision_models.md",
            "ionization_models.md",
            "anomalous_transport.md",
            "source_terms.md",
        ],
        "Internals" => "internals.md"
    ],
)

deploydocs(
    repo = "github.com/um-pepl/HallThruster.jl.git",
    devbranch = "main"
)