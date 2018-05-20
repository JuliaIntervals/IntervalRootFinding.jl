using IntervalArithmetic, StaticArrays, IntervalRootFinding

@testset "Linear Equations" begin

    A = [[2..3 0..1;1..2 2..3], ]
    b = [[0..120, 60..240], ]
    x = [[-120..90, -60..240], ]


    for i in 1:length(A)
        for precondition in (false, true)
            for f in (gauss_seidel_interval, gauss_seidel_contractor, gauss_elimination_interval)
                @test all(x[i] .⊆ f(A[i], b[i]))
            end
        end
    end
end
