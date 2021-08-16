using Documenter, HallThruster

push!(LOAD_PATH,"../src/")

makedocs(
    sitename = "HallThruster.jl",
    pages=[
        "Home" => "index.md",
        "Internals" => "internals.md"
    ],
)

deploydocs(
    repo = "github.com/archermarx/HallThruster.jl.git",
    devbranch = "main"
)