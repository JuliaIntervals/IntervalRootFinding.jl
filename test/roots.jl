
using IntervalArithmetic, IntervalRootFinding, StaticArrays, ForwardDiff
using Test

function all_unique(rts)
    all(root_status.(rts) .== :unique)
end

function test_newtonlike(f, deriv, X, method, nsol, tol=1e-3)
    rts = roots(f, X, method)
    @test length(rts) == nsol
    @test all_unique(rts)
    @test rts == roots(f, deriv, X, method)
end

newtonlike_methods = [Newton, Krawczyk]

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

    for method in newtonlike_methods
        test_newtonlike(sin, cos, -5..5, method, 3)
    end

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

    for method in newtonlike_methods
        deriv = xx -> ForwardDiff.jacobian(f, xx)
        test_newtonlike(f, deriv, X, method, 2)
    end

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

    for method in newtonlike_methods
        rts = roots(g, XX, method)
        @test length(rts) == 4
        @test all_unique(rts)
        deriv = xx -> ForwardDiff.jacobian(g, xx)
        @test rts == roots(g, deriv, XX, method)
    end
end

@testset "Stationary points" begin
    f(xx) = ( (x, y) = xx; sin(x) * sin(y) )
    gradf = ∇(f)
    XX = IntervalBox(-5..6, 2)
    tol = 1e-5

    for method in newtonlike_methods
        deriv = xx -> ForwardDiff.jacobian(gradf, xx)
        test_newtonlike(gradf, deriv, XX, method, 25, tol)
    end
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

    for method in newtonlike_methods
        deriv = z -> 3*z^2
        test_newtonlike(f, deriv, Xc, method, 3)
    end
end

@testset "RootSearch interface" begin
    contractor = Newton(sin, cos)
    search = RootSearch(-10..10, contractor, 1e-3)
    state, _ = iterate(search)

    # check that original and copy are independent
    state_copy = copy(state)
    pop!(state_copy.working) # mutate copy
    @test length(state.working) != length(state_copy.working)

    # cover optional iterator methods
    @test eltype(search) != Any
    # @test_nowarn iteratorsize(search)
end
