using IntervalArithmetic, IntervalRootFinding
using ForwardDiff
using Base.Test


const Dual = ForwardDiff.Dual

@testset "Promotion between Dual and Interval" begin
    @test promote(ForwardDiff.Dual(2, 1), Interval(1, 2)) ==
        (Dual(2..2, 1..1), Dual(1..2, 0..0))

    @test promote(Interval(2, 3), Dual(2, 1)) == (Dual(Interval(2, 3), Interval(0)),
            Dual(Interval(2.0), Interval(1.0)))

end
