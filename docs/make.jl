using Documenter
using IntervalArithmetic, IntervalRootFinding

makedocs(
    modules = [IntervalRootFinding],
    doctest = true,
    format = Documenter.HTML(),
    sitename = "IntervalRootFinding.jl",

    pages = Any[
        "Home" => "index.md",
        "`roots` interface" => "roots.md",
        "Internals" => "internals.md",
        "API" => "api.md"
        ]
    )

# deploydocs(
#     repo = "github.com/JuliaIntervals/IntervalRootFinding.jl.git",
#     target = "build",
#     deps = nothing,
#     make = nothing,
#     julia = "release",
#     osname = "linux"
# )
