using Documenter
using IntervalArithmetic, IntervalRootFinding

makedocs(
    modules = [IntervalRootFinding],
    doctest = true,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "IntervalRootFinding.jl",
    authors = "David P. Sanders and Luis Benet",
    pages = Any[
        "Home" => "index.md",
        "`roots` interface" => "roots.md",
        "Decorations and guarantees" => "decorations.md",
        "Internals" => "internals.md",
        "Bibliography" => "biblio.md",
        "API" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaIntervals/IntervalRootFinding.jl.git",
    target = "build",
    deps = nothing,
    make = nothing
)
