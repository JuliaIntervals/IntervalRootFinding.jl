using Test

import ForwardDiff

@testset "Stationary points" begin
    # Rastrigin function
    function rastrigin(x, y, A)
        return 2A + x^2 - A*cos(2π*x) + y^2 - A*cos(2π*y)
    end

    f(xy) = rastrigin(xy..., 10)

    ∇f = X -> ForwardDiff.gradient(f, X)

    L = 5.0
    X = interval(-L, (L+1))
    XX = [A, A]

    rts = IntervalRootFinding.roots(∇f, XX ; contractor = Newton, abstol = 1e-5)
    rts2 = IntervalRootFinding.roots(∇f, XX ; contractor = Krawczyk, abstol = 1e-5)

    @test length(rts) == length(rts2) == 529
    @test all(isunique, rts)
    @test all(isunique, rts2)
end