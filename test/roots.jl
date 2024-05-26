
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

newtonlike_methods = [Newton, Krawczyk]

@testset "1D roots" begin
    # Default
    rts = roots(sin, interval(-5, 5))
    @test length(rts) == 3
    @test all_unique(rts)

    # Bisection
    rts = roots(sin, interval(-5, 6) ; contractor = Bisection, abstol = 1e-3)
    @test length(rts) == 3

    # Refinement
    for rt in rts
        refined = roots(sin, rt ; contractor = Newton)
        @test length(refined) == 1
        @test isunique(only(refined))
    end

    for method in newtonlike_methods
        test_newtonlike(sin, cos, interval(-5, 5), method, 3)
    end

    # Infinite interval
    rts = roots(x -> x^2 - 2, interval(-Inf, Inf)) 
    @test length(rts) == 2
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

@testset "Stationary points" begin
    f(xx) = ( (x, y) = xx; sin(x) * sin(y) )
    gradf = xx -> ForwardDiff.gradient(f, xx)

    XX = [interval(-5, 6), interval(-5, 6)]
    tol = 1e-5

    for method in newtonlike_methods
        deriv = xx -> ForwardDiff.jacobian(gradf, xx)
        test_newtonlike(gradf, deriv, XX, method, 25, tol)
    end
end

@testset "Complex roots" begin
    @test_broken false
    return
    X = interval(-5, 5)
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
        test_newtonlike(f, deriv, Xc, method, 3, 1e-10)
    end
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