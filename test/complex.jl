using IntervalArithmetic, IntervalRootFinding
using Base.Test

@testset "Complex bisection" begin
    f(z) = z^4 + z - 2

    L = 10
    g, roots = complex_bisection(f, -L..L, -L..L)

    @test length(roots) == 8
end
