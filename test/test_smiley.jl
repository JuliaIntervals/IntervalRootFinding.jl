include("../examples/smiley_examples.jl")

using Test
using IntervalArithmetic, IntervalRootFinding
using .SmileyExample22, .SmileyExample52, .SmileyExample54, .SmileyExample55

function test_all_unique(xs)
    for x in xs
        @test x.status == :unique
    end
    return nothing
end

const tol = 1e-6
const method = Newton # NOTE: Bisection method performs badly in all examples

@info("Testing method $(method)")

@testset "$(SmileyExample22.title)" begin
    roots_found = roots(SmileyExample22.f, SmileyExample22.region, method, tol)
    @test length(roots_found) == 8
    test_all_unique(roots_found)
    # no reference data for roots given
end

for example in (SmileyExample52, SmileyExample54) #, SmileyExample55)
    @testset "$(example.title)" begin
        roots_found = roots(example.f, example.region, method, tol)
        @test length(roots_found) == length(example.known_roots)
        test_all_unique(roots_found)
        for rf in roots_found
            # check there is exactly one known root for each found root
            @test sum(!isempty(rk ∩ rf.interval)
                    for rk in example.known_roots) == 1
        end
    end
end
