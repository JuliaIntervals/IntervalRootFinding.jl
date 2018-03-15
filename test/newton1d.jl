
using IntervalArithmetic, IntervalRootFinding
using ForwardDiff
using Base.Test

const D = IntervalRootFinding.derivative

setprecision(Interval, Float64)
float_pi = @interval(pi)

setprecision(Interval, 10000)
big_pi = @interval(pi)
# Using precision "only" 256 leads to overestimation of the true roots for `cos`
# i.e the Newton method gives more accurate results!

half_pi = big_pi / 2
three_halves_pi = 3*big_pi/2

@testset "Testing newton1d" begin

    rts =  newton1d(sin, cos, -5..5)
    @test length(rts) == 3
    @test (-π.. -π) == rts[1].interval && :unique == rts[1].status
    @test (0.. 0) == rts[2].interval && :unique == rts[2].status
    @test (π.. π) == rts[3].interval && :unique == rts[3].status

    f(x) = e^(x^2) - cos(x)
    f′(x) = 2*x*e^(x^2) + sin(x)
    rts = newton1d(f, f′, -∞..∞)
    @test length(rts) == 1
    @test (0.. 0) == rts[1].interval && :unique == rts[1].status

    f(x) = x^4 - 10x^3 + 35x^2 - 50x + 24
    f′(x) = 4x^3 - 30x^2 + 70x - 50
    rts = newton1d(f, f′, -10..10)
    @test length(rts) == 4
    @test 1 in rts[1].interval && :unique == rts[1].status
    @test 2 in rts[2].interval && :unique == rts[2].status
    @test 3 in rts[3].interval && :unique == rts[3].status
    @test 4 in rts[4].interval && :unique == rts[4].status
end
