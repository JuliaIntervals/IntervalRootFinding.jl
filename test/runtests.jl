using IntervalRootFinding
using IntervalArithmetic.Symbols
using BranchAndPrune
using Test

newtonlike_methods = [Newton, Krawczyk]

include("roots.jl")
include("stationary_points.jl")
include("svectors.jl")
include("test_smiley.jl")
include("linear_eq.jl")

# include("newton1d.jl")
# include("quadratic.jl")
# include("slopes.jl")
