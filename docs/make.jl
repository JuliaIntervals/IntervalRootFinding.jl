using Documenter
using IntervalArithmetic, IntervalRootFinding

makedocs(
    modules = [IntervalRootFinding],
    doctest = true,
    format = :html,
    authors = "David P. Sanders",
    sitename = "IntervalRootFinding.jl",

    pages = Any[
        "Home" => "index.md",
        "API" => "api.md"
    ]
    )

# deploydocs(
#     deps = Deps.pip("pygments", "mkdocs", "mkdocs-cinder", "python-markdown-math"),
#     repo   = "github.com/dpsanders/IntervalConstraintProgramming.jl.git",
#     julia = "0.5",
#     osname = "linux"
# )

deploydocs(
    repo = "github.com/JuliaIntervals/IntervalRootFinding.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
    julia = "release",
    osname = "linux"
)
