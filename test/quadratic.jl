using IntervalArithmetic, IntervalRootFinding

@testset "Quadratic Interval Equations" begin
    a = [interval(1, 5), interval(-1, 1), interval(1, 2), ]
    b = [interval(2, 4), interval(-2, 2), interval(1, 3), ]
    c = [interval(0, 2), interval(-1, 1), interval(2, 6), ]
    x = [[interval(-4, -2), interval(-0, -0)], [interval(-∞, -1), interval(-1, -1), interval(-1, ∞)], [interval(-2, -1)], ]
    l = [2, 3, 1, ]
    for i in 1:length(a)
            roots = quadratic_interval(a[i], b[i], c[i])
            @test length(roots) == l[i]
            @test all(roots == x[i])
    end
end
