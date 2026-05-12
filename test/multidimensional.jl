@testset "2D roots" begin
    f(x, y) = [x^2 + y^2 - 1, y - 2x]
    f(X) = f(X...)
    X = [interval(-6, 6), interval(-6, 6)]

    # Bisection
    rts = roots(f, X ; contractor = Bisection, abstol = 1e-3)
    exact_sol = [sqrt(1/5), 2sqrt(1/5)]
    @test in_solution_set(exact_sol, root_region.(rts))
    @test in_solution_set(-exact_sol, root_region.(rts))

    for contractor in newtonlike_methods
        deriv = xx -> ForwardDiff.jacobian(f, xx)
        test_newtonlike(f, deriv, X, contractor, 2)
        rts = roots(f, X ; contractor)
        @test in_solution_set(exact_sol, root_region.(rts))
        @test in_solution_set(-exact_sol, root_region.(rts))
    end

    # Infinite interval
    X = [interval(-Inf, Inf), interval(-Inf, Inf)]
    @suppress rts = roots(f, X ; contractor = Newton)
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
    dg = xx -> ForwardDiff.jacobian(g, xx)

    X = interval(-5, 5)
    XX = [X, X, X]

    for contractor in newtonlike_methods
        rts = roots(g, XX ; contractor)
        @test length(rts) == 4
        @test all(isunique, rts)
        @test all(rts .== roots(g, XX ; contractor, derivative = dg))
    end
end

@testset "10D roots" begin
    include("../examples/10_dimensional.jl")

    for contractor in newtonlike_methods
        rts = roots(f10d, X10d_large ; contractor, max_iteration = 50_000)
        @test any(X -> all(in_interval.(sol10d, X.region)), rts)

        rts = roots(f10d, X10d_close ; contractor, max_iteration = 1_000_000)
        @test length(rts) == 1
        @test isunique(rts[1])
    end
end

@testset "Dimension mismatch" begin
    f21(xy) = [xy[1]^2 - 2]
    f23(xy) = [xy[1]^2 - 2, xy[2]^2 - 3, xy[1] + xy[2]]

    X = [interval(0, 5), interval(0, 5)]

    for contractor in newtonlike_methods
        @test_throws DimensionMismatch roots(f21, X ; contractor)
        @test_throws DimensionMismatch roots(f23, X ; contractor)
    end
end
