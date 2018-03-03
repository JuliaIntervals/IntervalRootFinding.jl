
using IntervalArithmetic, IntervalRootFinding, StaticArrays
using Base.Test

@testset "1D roots" begin
    rts = roots(sin, -5..5)
    @test length(rts) == 4
    @test length(find(x->x==:unique, [root.status for root in rts])) == 2

    rts = roots(sin, -5..6, Bisection)
    @test length(rts) == 3
    rts = roots(sin, rts, Newton)
    @test all(root.status == :unique for root in rts)
end


@testset "2D roots" begin
    f(x, y) = SVector(x^2 + y^2 - 1, y - 2x)
    f(X) = f(X...)

    rts = roots(f, (-6..6) × (-6..6), Bisection, 1e-3)
    @test length(rts) == 4

    rts = roots(f, rts, Newton)
    @test length(rts) == 2
end

@testset "Complex roots" begin
    x = -5..6
    Xc = Complex(x, x)
    f(z) = z^3 - 1

    rts = roots(f, Xc, Bisection, 1e-3)
    @test length(rts) == 7
    rts = roots(f, rts, Newton)
    @test length(rts) == 3
    rts = roots(f, Xc)
    @test length(rts) == 3
end

# From R docs:
# https://www.rdocumentation.org/packages/pracma/versions/1.9.9/topics/broyden

@testset "3D roots" begin
    function g(x)
        (x1, x2, x3) = x
        SVector(    x1^2 + x2^2 + x3^2 - 1,
                    x1^2 + x3^2 - 0.25,
                    x1^2 + x2^2 - 4x3
                )
    end

    X = (-5..5)
    rts = roots(g, IntervalBox(X, 3))
    @test length(rts) == 4
end
