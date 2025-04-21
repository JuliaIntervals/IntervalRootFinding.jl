include("../examples/smiley_examples.jl")
using .SmileyExample22, .SmileyExample52, .SmileyExample54, .SmileyExample55

using BenchmarkTools
using ForwardDiff
using IntervalArithmetic
using IntervalRootFinding
using StaticArrays

import Random

const SUITE = BenchmarkGroup()

Random.seed!(0)  # Seed the RNG to get consistent results
tol = 1e-10

include("dietmar_ratz_functions.jl")

S = SUITE["Smiley"] = BenchmarkGroup()
for example in (SmileyExample22, SmileyExample52, SmileyExample54) #, SmileyExample55)
    s = S[example.title] = BenchmarkGroup()
    for contractor in (Newton, Krawczyk)
        s[string(contractor)] = @benchmarkable roots($(example.f), $(example.region) ; contractor = $contractor, abstol = $tol)
    end
end


S = SUITE["Rastigrin stationary points"] = BenchmarkGroup()

# Rastrigin function:
const A = 10

f(x, y) = 2A + x^2 - A*cos(2π*x) + y^2 - A*cos(2π*y)
f(X) = f(X...)

∇f = X -> ForwardDiff.gradient(f, X)

L = 5.0
X = SVector(interval(-L, (L+1)), interval(-L, (L+1)))

for contractor in (Newton, Krawczyk)
    S[string(contractor)] = @benchmarkable roots($(∇f), $X ; contractor = $contractor, abstol = 1e-5)
end


S = SUITE["Linear equations"] = BenchmarkGroup()

sizes = (2, 5, 10)

for n in sizes
    s = S["n = $n"] = BenchmarkGroup()
    M = interval.(randn(n, n))
    b = interval.(randn(n))

    # s["Gauss seidel"] = @benchmarkable gauss_seidel_interval($M, $b)
    # s["Gauss seidel contractor"] = @benchmarkable gauss_seidel_contractor($M, $b)
    # s["Gauss elimination"] = @benchmarkable gauss_elimination_interval($M, $b)
end


S = SUITE["Dietmar-Ratz"] = BenchmarkGroup()
X = interval(0.75, 1.75)

for (k, dr) in enumerate(dr_functions)
    s = S["Dietmar-Ratz $k"] = BenchmarkGroup()

    if k != 8  # dr8 is excluded as it has too many roots
        for contractor in (Newton, Krawczyk)
            s[string(contractor)] = @benchmarkable roots($dr, $X ; contractor = $contractor, abstol = $tol)
        end
    end
end
