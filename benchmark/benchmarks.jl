include("../examples/smiley_examples.jl")
using .SmileyExample22, .SmileyExample52, .SmileyExample54, .SmileyExample55

using BenchmarkTools
using IntervalRootFinding

const SUITE = BenchmarkGroup()

tol = 1e-15

SUITE["Smiley"] = BenchmarkGroup()
for method in (Newton, Krawczyk)
    SUITE["Smiley"][string(method)] = BenchmarkGroup()
    for example in (SmileyExample22, SmileyExample52, SmileyExample54) #, SmileyExample55)
        SUITE["Smiley"][string(method)][example.title] = @benchmarkable roots($(example.f), $(example.region), $method, $tol)
    end
end
