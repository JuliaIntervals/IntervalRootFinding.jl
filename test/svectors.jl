@testset "Static vectors" begin
    f(xx) = [xx[1]^2 - 1, xx[2]^2 - 2]
    g(xx) = SVector(xx[1]^2 - 1, xx[2]^2 - 2)

    X = [interval(-5, 5), interval(-5, 5)]
    S = SVector(interval(-5, 5), interval(-5, 5))

    @testset for contractor in newtonlike_methods
        rts = roots(f, X ; contractor)
        rts2 = roots(g, S ; contractor)

        for (rt, rt2) in zip(rts, rts2)
            @test isequal_interval(
                root_region(rt),
                root_region(rt2)
            )
        end
    end
end