using IntervalArithmetic, StaticArrays, IntervalRootFinding
using Test

function rand_vec(n::Int)
    a = randn(n)
    A = interval.(a)
    mA = MVector{n}(A)
    sA = SVector{n}(A)
    return A, mA, sA
end

function rand_mat(n::Int)
    a = randn(n, n)
    A = interval.(a)
    mA = MMatrix{n, n}(A)
    sA = SMatrix{n, n}(A)
    return A, mA, sA
end

@testset "Linear Equations" begin

    As = [[2..3 0..1; 1..2 2..3], ]
    bs = [[0..120, 60..240], ]
    xs = [[-120..90, -60..240], ]

    for i in 1:10
        rand_A = rand_mat(i)[1]
        rand_x = rand_vec(i)[1]
        rand_b = rand_A * rand_x
        push!(As, rand_A)
        push!(bs, rand_b)
        push!(xs, rand_x)
    end

    n = length(As)

    for solver in (gauss_seidel_interval, gauss_seidel_contractor, gauss_elimination_interval, \)
        @testset "Solver $solver" begin
            if solver in (gauss_seidel_interval, gauss_seidel_contractor)
                @test_broken false
                continue
            end

            for i in 1:n
                soln = solver(As[i], bs[i])
                @test all(issubset_interval.(xs[i], soln))
            end
        end
    end
end
