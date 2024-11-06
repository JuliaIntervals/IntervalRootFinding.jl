using Test

import ForwardDiff

@testset "Stationary points" begin
    f(xx) = ( (x, y) = xx; sin(x) * sin(y) )
    gradf = xx -> ForwardDiff.gradient(f, xx)

    XX = [interval(-5, 6), interval(-5, 6)]
    tol = 1e-5

    for method in newtonlike_methods
        deriv = xx -> ForwardDiff.jacobian(gradf, xx)
        test_newtonlike(gradf, deriv, XX, method, 25, tol)
    end

    # Rastrigin function
    function rastrigin(x, y, A)
        return 2A + x^2 - A*cos(2π*x) + y^2 - A*cos(2π*y)
    end

    f(xy) = rastrigin(xy..., 10)

    ∇f = X -> ForwardDiff.gradient(f, X)

    L = 5.0
    X = interval(-L, (L+1))
    XX = [X, X]

    rts = IntervalRootFinding.roots(∇f, XX ; contractor = Newton, abstol = 1e-5)
    rts2 = IntervalRootFinding.roots(∇f, XX ; contractor = Krawczyk, abstol = 1e-5)

    @test length(rts) == length(rts2) == 529
    @test all(isunique, rts)
    @test all(isunique, rts2)
end