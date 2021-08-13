using Documenter, HallThruster

push!(LOAD_PATH,"../src/")
makedocs(sitename="HallThruster.jl Documentation")

deploydocs(
    repo = "github.com/archermarx/HallThruster.jl.git",
)