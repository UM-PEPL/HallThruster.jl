using Documenter, HallThruster

push!(LOAD_PATH,"../src/")

makedocs(
    modules=[HallThruster],
    authors="Thomas Marks <marksta@umich.edu> and contributors",
    repo="https://github.com/archermarx/HallThruster.jl/blob/{commit}{path}#L{line}",
    sitename="HallThruster.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://archermarx.github.io/HallThruster.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Internals" => "internals.md"
    ],
)

deploydocs(
    repo = "github.com/archermarx/HallThruster.jl.git",
    devbranch = "main"
)