using ValidatedNumerics, ValidatedNumerics.RootFinding
using Base.Test

@testset "Bisection tests" begin
    @testset "Single variable" begin
        f(x) = sin(x) + cos(x/2) + x/5

        roots = bisection(f, -∞..∞, tolerance=1e-3)
        roots2 = [root.interval for root in roots]

        @test roots2 ==  [  Interval(-0.8408203124999999, -0.8398437499999999),
                            Interval(3.6425781249999996, 3.6435546874999996),
                            Interval(6.062499999999999, 6.063476562499999)
                         ]

    end

    @testset "Two variables" begin
        @intervalbox g(x, y) = (x^2 + y^2 - 1, y - x)

        X = (-∞..∞) × (-∞..∞)
        roots = bisection(g, X, tolerance=1e-3)
        roots2 = [root.interval for root in roots]
        
        @test roots2 == [IntervalBox(Interval(-0.7080078124999999, -0.7070312499999999), Interval(-0.7080078124999999, -0.7070312499999999)),
 IntervalBox(Interval(-0.7080078124999999, -0.7070312499999999), Interval(-0.7070312499999999, -0.7060546874999999)),
 IntervalBox(Interval(-0.7070312499999999, -0.7060546874999999), Interval(-0.7080078124999999, -0.7070312499999999)),
 IntervalBox(Interval(0.7060546874999999, 0.7070312499999999), Interval(0.7070312499999999, 0.7080078124999999)),
 IntervalBox(Interval(0.7070312499999999, 0.7080078124999999), Interval(0.7060546874999999, 0.7070312499999999)),
 IntervalBox(Interval(0.7070312499999999, 0.7080078124999999), Interval(0.7070312499999999, 0.7080078124999999))]

    end

end
