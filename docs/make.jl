using Documenter, HallThruster

push!(LOAD_PATH,"../src/")

makedocs(
    pages=[
        "Home" => "index.md",
        "Internals" => "internals.md"
    ],
)

deploydocs(
    repo = "github.com/archermarx/HallThruster.jl.git",
    devbranch = "main"
)