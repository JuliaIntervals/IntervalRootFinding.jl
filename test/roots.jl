
using IntervalArithmetic, IntervalRootFinding, StaticArrays
using Base.Test

@testset "1D roots" begin
    rts = roots(sin, -5..5)
    @test length(rts) == 3
    @test length(find(x->x==:unique, [root.status for root in rts])) == 3

    rts = roots(sin, -5..6, Bisection)
    @test length(rts) == 3
    rts = roots(sin, rts, Newton)
    @test all(root.status == :unique for root in rts)

    rts = roots(sin, -5..5, Newton)
    @test rts == roots(sin, -5..5, Newton; deriv = cos)
end


@testset "2D roots" begin
    f(x, y) = SVector(x^2 + y^2 - 1, y - 2x)
    f(X) = f(X...)
    X = (-6..6) × (-6..6)

    rts = roots(f, X, Bisection, 1e-3)
    @test length(rts) == 4

    rts = roots(f, rts, Newton)
    @test length(rts) == 2

    rts = roots(f, X, Newton)
    @test rts == roots(f, X, Newton; deriv = xx -> ForwardDiff.jacobian(f, xx))
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
    XX = IntervalBox(X, 3)
    rts = roots(g, XX, Newton)
    @test length(rts) == 4
    @test rts == roots(g, XX, Newton; deriv = xx -> ForwardDiff.jacobian(g, xx))
end

@testset "Stationary points" begin
    f(xx) = ( (x, y) = xx; sin(x) * sin(y) )
    gradf = ∇(f)
    XX = IntervalBox(-5..6, 2)
    tol = 1e-5

    rts = roots(gradf, XX, Newton, tol)
    @test length(rts) == 25
    @test rts == roots(gradf, XX, Newton, tol; deriv = xx -> ForwardDiff.jacobian(gradf, xx))
end

@testset "Complex roots" begin
    X = -5..5
    Xc = Complex(X, X)
    f(z) = z^3 - 1

    rts = roots(f, Xc, Bisection, 1e-3)
    @test length(rts) == 7
    rts = roots(f, rts, Newton)
    @test length(rts) == 3
    rts = roots(f, Xc)
    @test length(rts) == 3

    rts = roots(f, Xc, Newton)
    rts2 = roots(f, Xc, Newton; deriv = z -> 3*z^2)

    intervals = [reim(rt.interval) for rt in rts]
    intervals2 = [reim(rt.interval) for rt in rts2]

    d = []
    for (I, I2) in zip(sort(intervals), sort(intervals2))
        append!(d, abs.(I .- I2))
    end
    @test all(d .< 1e-15)
end
