
using IntervalArithmetic, IntervalRootFinding, StaticArrays
using Base.Test

function all_unique(rts)
    all(root_status.(rts) .== :unique)
end

@testset "1D roots" begin
    # Default
    rts = roots(sin, -5..5)
    @test length(rts) == 3
    @test all_unique(rts)

    # Bisection
    rts = roots(sin, -5..6, Bisection)
    @test length(rts) == 3

    # Refinement
    rts = roots(sin, rts, Newton)
    @test all_unique(rts)

    # Newton
    rts = roots(sin, -5..5, Newton)
    @test length(rts) == 3
    @test all_unique(rts)
    @test rts == roots(sin, -5..5, Newton; deriv = cos)

    # Krawczyz
    rts = roots(sin, -5..5, Krawczyk)
    @test length(rts) == 3
    @test all_unique(rts)
    @test rts == roots(sin, -5..5, Krawczyk; deriv = cos)

    # Infinite interval
    rts = roots(x -> x^2 - 2, -∞..∞)
    @test length(rts) == 2
end


@testset "2D roots" begin
    f(x, y) = SVector(x^2 + y^2 - 1, y - 2x)
    f(X) = f(X...)
    X = (-6..6) × (-6..6)

    # Bisection
    rts = roots(f, X, Bisection, 1e-3)
    @test length(rts) == 4

    # Newton
    rts = roots(f, rts, Newton)
    @test length(rts) == 2
    @test all_unique(rts)

    rts = roots(f, X, Newton)
    @test rts == roots(f, X, Newton; deriv = xx -> ForwardDiff.jacobian(f, xx))

    # Krawczyk
    rts = roots(f, X, Krawczyk)
    @test length(rts) == 2
    @test all_unique(rts)
    @test rts == roots(f, X, Krawczyk; deriv = xx -> ForwardDiff.jacobian(f, xx))

    # Infinite interval
    X = IntervalBox(-∞..∞, 2)
    rts = roots(f, X, Newton)
    @test length(rts) == 2
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

    # Newton
    rts = roots(g, XX, Newton)
    @test length(rts) == 4
    @test all_unique(rts)
    @test rts == roots(g, XX, Newton; deriv = xx -> ForwardDiff.jacobian(g, xx))

    # Krawczyk
    rts = roots(g, XX, Krawczyk)
    @test length(rts) == 4
    @test all_unique(rts)
    @test rts == roots(g, XX, Krawczyk; deriv = xx -> ForwardDiff.jacobian(g, xx))
end

@testset "Stationary points" begin
    f(xx) = ( (x, y) = xx; sin(x) * sin(y) )
    gradf = ∇(f)
    XX = IntervalBox(-5..6, 2)
    tol = 1e-5

    # Newton
    rts = roots(gradf, XX, Newton, tol)
    @test length(rts) == 25
    @test all_unique(rts)
    @test rts == roots(gradf, XX, Newton, tol; deriv = xx -> ForwardDiff.jacobian(gradf, xx))

    # Krawczyk
    rts = roots(gradf, XX, Krawczyk, tol)
    @test length(rts) == 25
    @test all_unique(rts)
    @test rts == roots(gradf, XX, Krawczyk, tol; deriv = xx -> ForwardDiff.jacobian(gradf, xx))
end

@testset "Complex roots" begin
    X = -5..5
    Xc = Complex(X, X)
    f(z) = z^3 - 1

    # Default
    rts = roots(f, Xc)
    @test length(rts) == 3

    # Bisection
    rts = roots(f, Xc, Bisection, 1e-3)
    @test length(rts) == 7

    # Newton
    rts = roots(f, rts, Newton)
    @test length(rts) == 3
    @test all_unique(rts)

    rts = roots(f, Xc, Newton)
    @test rts == roots(f, Xc, Newton; deriv = z -> 3*z^2)

    # Krawczyk
    rts = roots(f, Xc, Krawczyk)
    @test length(rts) == 3
    @test all_unique(rts)
    @test rts == roots(f, Xc, Krawczyk; deriv = z -> 3*z^2)
end
