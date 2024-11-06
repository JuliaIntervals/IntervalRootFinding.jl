
using IntervalArithmetic, IntervalRootFinding, StaticArrays, ForwardDiff
using Test

function all_unique(rts)
    all(root_status.(rts) .== :unique)
end

function roots_dist(rt1::Root{<:Interval}, rt2::Root{<:Interval})
    d = dist(root_region(rt1), root_region(rt2))
    return sum(d)
end

function roots_dist(rt1::Root{<:AbstractVector}, rt2::Root{<:AbstractVector})
    return sum(dist.(root_region(rt1), root_region(rt2)))
end

function roots_dist(rt1::Root{Complex{<:Interval}}, rt2::Root{Complex{<:Interval}})
    dreal = dist(real(root_region(rt1)), real(root_region(rt2)))
    dimag = dist(imag(root_region(rt1)), imag(root_region(rt2)))

    return sum(dreal) + sum(dimag)
end

function test_newtonlike(f, derivative, X, contractor, nsol, tol=1e-10)
    rts = roots(f, X ; contractor)
    @test length(rts) == nsol
    @test all_unique(rts)
    @test sum(roots_dist.(rts, roots(f, X ; contractor, derivative))) < tol
end

@testset "1D roots" begin
    # Default
    rts = roots(sin, interval(-5, 5))
    @test length(rts) == 3
    @test all_unique(rts)

    # Bisection
    rts = roots(sin, interval(-5, 6) ; contractor = Bisection, abstol = 1e-3)
    @test length(rts) == 3

    @testset for contractor in newtonlike_methods
        # Refinement
        rts = roots(sin, interval(-5, 6) ; contractor = Bisection, abstol = 1e-3)
        for rt in rts
            refined = roots(sin, rt ; contractor)
            @test length(refined) == 1
            @test isunique(only(refined))
        end

        test_newtonlike(sin, cos, interval(-5, 5), contractor, 3)

        # Infinite interval
        rts = roots(x -> x^2 - 2, interval(-Inf, Inf) ; contractor)
        @test length(rts) == 2

        # abs
        g(p) = abs(1 / ( (1+p)^30 ) * 10_000 - 100)
        g2(p) = (1 / ( (1+p)^30 ) * 10_000 - 100)^2
        X = interval(-1000, 1000)
        rts = roots(g, X ; contractor)
        rts2 = roots(g2, X ; contractor)
        rr = intersect_interval.(root_region.(rts), root_region.(rts2))
        @test !any(isempty(rr))

        # Hard problem with a known zero at zero
        h(E) = tan(√(2E)) + √(E / (3π^2 - E))
        rts = roots(h, interval(-1e100, 1e100) ; contractor)
        @test any(rts) do rt
            in_interval(0, root_region(rt))
        end

        # Singularity
        s(x) = x + zero(x)/x
        rts = roots(s, interval(-1, 1) ; contractor)
        @test length(rts) == 1
        @test root_status(only(rts)) == :unknown
        @test in_interval(0, root_region(only(rts)))
    end
end


function in_solution_set(point, solution_intervals)
    return any(map(Y -> in_region(point, Y), solution_intervals))
end

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
    rts = roots(f, X ; contractor = Newton)
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
        @test all_unique(rts)
        @test all(rts .== roots(g, XX ; contractor, derivative = dg))
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

@testset "Out of domain" begin
    for contractor in newtonlike_methods
        @test length(roots(log, interval(-100, 2) ; contractor)) == 1
        @test length(roots(log, interval(-100, -1) ; contractor)) == 0
    end
end

@testset "Infinite domain" begin
    for contractor in newtonlike_methods
        rts = roots(x -> x^2 - 2, interval(-Inf, Inf) ; contractor)
        @test length(filter(isunique, rts)) == 2
    end
end

@testset "NaN return value" begin
    f(xx) = ( (x, y) = xx; [log(y/x) + 3x, y - 2x] )
    X = [interval(-100, 100), interval(-100, 100)]
    for contractor in newtonlike_methods
        rts = roots(f, X ; contractor)
        @test length(filter(isunique, rts)) == 1
        @test length(filter(x -> all(in_interval.(0, x)), root_region.(rts))) == 1
    end
end

@testset "Root at infinity" begin
    for contractor in newtonlike_methods
        pb = RootProblem(x -> 1/x, interval(1, Inf) ; contractor)

        state = nothing
        for (k, s) in enumerate(pb)
            state = s
            
            k > 1000 && break
        end

        rts = [leaf.region for leaf in BranchAndPrune.Leaves(state.tree)]
        for rt in rts
            # Some roots have an empty trivial interval for unknown reason,
            # which is unhelpful but not incorrect
            if decoration(rt.region) != trv
                @test sup(rt.region) == Inf
            end
        end
    end
end

@testset "Complex roots" begin
    X = interval(-5, 5)
    Xc = Complex(X, X)
    f(z) = z^3 - 1
    df(z) = 3z^2

    # Default
    for contractor in newtonlike_methods
        rts = roots(f, Xc ; derivative = df, contractor, abstol = 1e-3)
        @test count(isunique, rts) == 3
    end

    # Bisection
    rts = roots(f, Xc ; contractor = Bisection, abstol = 1e-3)
end

@testset "@exact" begin
    @testset "1D" begin
        @exact f(x) = x^2 - 2.25
        @exact g(x) = sin(1.57x)

        f0(x) = x^2 - 2.25
        g0(x) = sin(1.57x)

        X = interval(-2.2, 2.2)

        @testset for method in newtonlike_methods
            for (func, func0, nsol) in [(f, f0, 2), (g, g0, 3)]
                rts = roots(func, X ; contractor = method)
                @test length(rts) == nsol

                rts0 = roots(func0, X ; contractor = method)
                for (rt, rt0) in zip(rts, rts0)
                    @test isequal_interval(root_region(rt), root_region(rt0))
                end

                @test all(isguaranteed.(root_region.(rts)))
            end
        end
    end

    @testset "2D" begin
        @exact f(xx) = [xx[1]^2 - 1, xx[2]^2 - 5]
        X = [interval(-2, 2), interval(-10, 10)]
        @testset for method in newtonlike_methods
            rts = roots(f, X ; contractor = method)
            @test all(rts) do rt
                all(isguaranteed.(root_region(rt)))
            end
        end
    end
end