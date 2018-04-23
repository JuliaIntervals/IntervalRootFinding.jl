
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

    f(x) = exp(x^2) - cos(x)
    f′(x) = 2*x*exp(x^2) + sin(x)
    f1(x) = x^4 - 10x^3 + 35x^2 - 50x + 24
    f1′(x) = 4x^3 - 30x^2 + 70x - 50

    for autodiff in (false, true)
        if autodiff
            rts1 = newton1d(sin, -5..5)
            rts2 = newton1d(f, -∞..∞)
            rts3 = newton1d(f1, -10..10)

        else
            rts1 = newton1d(sin, cos, -5..5)
            rts2 = newton1d(f, f′, -∞..∞)
            rts3 = newton1d(f1, f1′, -10..10)
        end

        @test length(rts1) == 3
        L = [ -pi_interval(Float64), 0..0, pi_interval(Float64)]
        for i = 1:length(rts1)
            @test L[i] == rts1[i].interval && :unique == rts1[i].status
        end

        @test length(rts2) == 1
        @test (0..0) == rts2[1].interval && :unique == rts2[1].status

        @test length(rts3) == 4
        L = [1, 2, 3, 4]
        for i = 1:length(rts3)
            @test L[i] in rts3[i].interval && :unique == rts3[i].status
        end
    end
end
