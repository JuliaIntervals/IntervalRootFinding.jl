using BenchmarkTools
using ForwardDiff
using IntervalArithmetic
using IntervalRootFinding
import Random
using StaticArrays

const SUITE = BenchmarkGroup()

Random.seed!(0)  # Seed the RNG to get consistent results
tol = 1e-10

## Smiley and Chun examples
include("../examples/smiley_examples.jl")
using .SmileyExample22, .SmileyExample52, .SmileyExample54, .SmileyExample55

S = SUITE["Smiley"] = BenchmarkGroup()

for example in (SmileyExample22, SmileyExample52, SmileyExample54, SmileyExample55)
    s = S[example.title] = BenchmarkGroup()
    for contractor in (Newton, Krawczyk)
        s[string(contractor)] = @benchmarkable roots($(example.f), $(example.region) ; contractor = $contractor, abstol = $tol, infer_root_type = false)
    end
end


## Rastrigin function
S = SUITE["Rastigrin stationary points"] = BenchmarkGroup()

const A = 10

f(x, y) = 2A + x^2 - A*cos(2π*x) + y^2 - A*cos(2π*y)
f(X) = f(X...)

∇f = X -> ForwardDiff.gradient(f, X)

L = 5.0
X = SVector(interval(-L, (L+1)), interval(-L, (L+1)))

for contractor in (Newton, Krawczyk)
    S[string(contractor)] = @benchmarkable roots($(∇f), $X ; contractor = $contractor, abstol = 1e-5, infer_root_type = false)
end


## Dietmar-Ratz functions
include("dietmar_ratz_functions.jl")

S = SUITE["Dietmar-Ratz"] = BenchmarkGroup()

X = interval(0.75, 1.75)

for (k, dr) in enumerate(dr_functions)
    s = S["Dietmar-Ratz $k"] = BenchmarkGroup()

    if k != 8  # dr8 is excluded as it has too many roots
        for contractor in (Newton, Krawczyk)
            s[string(contractor)] = @benchmarkable roots($dr, $X ; contractor = $contractor, abstol = $tol, infer_root_type = false)
        end
    end
end

# 10 dimensional problem
include("../examples/10_dimensional.jl")

S = SUITE["10 dimensional"] = BenchmarkGroup()

for contractor in (Newton, Krawczyk)
    S[string(contractor)] = @benchmarkable roots($f10d, $X10d_close ; contractor = $contractor, abstol = $tol, max_iteration = 1_000_000)
end