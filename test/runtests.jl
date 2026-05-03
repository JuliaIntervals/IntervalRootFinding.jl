using BranchAndPrune
using IntervalRootFinding
using IntervalArithmetic.Symbols
using ForwardDiff
using StaticArrays
using Suppressor
using Test

newtonlike_methods = [Newton, Krawczyk]

# Compute the distance between set of roots
function roots_dist(rt1::Root{<:Interval}, rt2::Root{<:Interval})
    d = dist(root_region(rt1), root_region(rt2))
    return sum(d)
end

function roots_dist(rt1::Root{<:AbstractVector}, rt2::Root{<:AbstractVector})
    return sum(dist.(root_region(rt1), root_region(rt2)))
end

function roots_dist(rt1::Root{Complex{<:Interval}}, rt2::Root{Complex{<:Interval}})
    dreal = dist(real(root_region(rt1)), real(root_region(rt2)))
    dimag = dist(imag(root_region(rt1)), imag(root_region(rt2)))

    return sum(dreal) + sum(dimag)
end

# Test that a contractor that use derivatives
# Check that all return roots are unique, and that giving the derivative explictly
# makes no difference
function test_newtonlike(f, derivative, X, contractor, nsol, tol=1e-10)
    rts = roots(f, X ; contractor)
    @test length(rts) == nsol
    @test all(isunique, rts)
    @test sum(roots_dist.(rts, roots(f, X ; contractor, derivative))) < tol
end

# Check that a given point is in one of the given intervals
function in_solution_set(point, solution_intervals)
    return any(map(Y -> in_region(point, Y), solution_intervals))
end

include("roots.jl")
include("stationary_points.jl")

include("multidimensional.jl")
include("svectors.jl")

include("test_smiley.jl")
include("linear_eq.jl")
