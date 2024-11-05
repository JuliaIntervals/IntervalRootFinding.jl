using IntervalRootFinding
using IntervalArithmetic.Symbols
using Test

newtonlike_methods = [Newton, Krawczyk]

include("roots.jl")
include("svectors.jl")
include("test_smiley.jl")
include("linear_eq.jl")

# include("newton1d.jl")
# include("quadratic.jl")
# include("slopes.jl")
