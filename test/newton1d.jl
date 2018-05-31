
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

    f(x) = exp(x^2) - cos(x)    #double root
    f′(x) = 2*x*exp(x^2) + sin(x)
    f1(x) = x^4 - 10x^3 + 35x^2 - 50x + 24  #four unique roots
    f1′(x) = 4x^3 - 30x^2 + 70x - 50
    f2(x) = 4567x^2 - 9134x + 4567  #double root
    f2′(x) = 9134x - 9134
    f3(x) = (x^2 - 2)^2 #two double roots
    f3′(x) = 4x * (x^2 - 2)
    f4(x) = sin(x) - x  #triple root at 0
    f4′(x) = cos(x) - 1
    f5(x) = (x^2 - 1)^4 * (x^2 - 2)^4 #two quadruple roots
    f5′(x) = 8x * (-3 + 2x^2) * (2 - 3x^2 + x^4)^3
    for autodiff in (false, true)
        if autodiff
            rts1 = newton1d(sin, -5..5)
            rts2 = newton1d(f, -∞..∞)
            rts3 = newton1d(f1, -10..10)
            rts4 = newton1d(f2, -10..11)
            rts5 = newton1d(f3, -10..10)
            rts6 = newton1d(f4, -10..10)
            rts7 = newton1d(f5, -10..10)

        else
            rts1 = newton1d(sin, cos, -5..5)
            rts2 = newton1d(f, f′, -∞..∞)
            rts3 = newton1d(f1, f1′, -10..10)
            rts4 = newton1d(f2, f2′, -10..11)
            rts5 = newton1d(f3, f3′, -10..10)
            rts6 = newton1d(f4, f4′, -10..10)
            rts7 = newton1d(f5, f5′, -10..10, reltol=0)
        end

        @test length(rts1) == 3
        L = [ -pi_interval(Float64), 0..0, pi_interval(Float64)]
        for i = 1:length(rts1)
            @test L[i] in rts1[i].interval && :unique == rts1[i].status
        end

        @test length(rts2) == 1
        @test (0..0) == rts2[1].interval && :unknown == rts2[1].status

        @test length(rts3) == 4
        L = [1, 2, 3, 4]
        for i = 1:length(rts3)
            @test L[i] in rts3[i].interval && :unique == rts3[i].status
        end

        @test length(rts4) == 1
        @test 1 in rts4[1].interval && :unknown == rts4[1].status

        L1 = [-sqrt(2), sqrt(2)]
        for i = 1:length(rts5)
            @test L1[i] in rts5[i].interval && :unknown == rts5[i].status
        end

        @test length(rts6) == 1
        @test 0 in rts6[1].interval && :unknown == rts6[1].status

        @test length(rts7) == 4
        L = [-sqrt(2), -1, 1, sqrt(2)]
        for i = 1:length(rts7)
            @test L[i] in rts7[i].interval && :unknown == rts7[i].status
        end

    end
end

@testset "Testing slope newton1d" begin

    f(x) = exp(x^2) - cos(x)    #double root
    f1(x) = x^4 - 10x^3 + 35x^2 - 50x + 24  #four unique roots
    f2(x) = 4567x^2 - 9134x + 4567  #double root
    f3(x) = (x^2 - 2)^2 #two double roots
    f4(x) = sin(x) - x  #triple root at 0
    f5(x) = (x^2 - 1)^4 * (x^2 - 2)^4 #two quadruple roots

    rts1 = slope_newton1d(sin, -5..5)
    rts2 = slope_newton1d(f, -∞..∞)
    rts3 = slope_newton1d(f1, -10..10)
    rts4 = slope_newton1d(f3, -10..10)
    rts5 = slope_newton1d(f4, -10..10)

    @test length(rts1) == 3
    L = [ -pi_interval(Float64), 0..0, pi_interval(Float64)]
    for i = 1:length(rts1)
        @test L[i] in rts1[i].interval && :unique == rts1[i].status
    end

    @test length(rts2) == 1
    @test (0..0) == rts2[1].interval && :unknown == rts2[1].status

    @test length(rts3) == 4
    L = [1, 2, 3, 4]
    for i = 1:length(rts3)
        @test L[i] in rts3[i].interval && :unique == rts3[i].status
    end

    L1 = [-sqrt(2), sqrt(2)]
    for i = 1:length(rts4)
        @test L1[i] in rts4[i].interval && :unknown == rts4[i].status
    end

    @test length(rts5) == 1
    @test 0 in rts5[1].interval && :unknown == rts5[1].status

end
