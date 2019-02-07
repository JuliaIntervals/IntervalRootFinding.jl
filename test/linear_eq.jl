using IntervalArithmetic, StaticArrays, IntervalRootFinding
using Test

function rand_vec(n::Int)
    a = randn(n)
    A = Interval.(a)
    mA = MVector{n}(A)
    sA = SVector{n}(A)
    return A, mA, sA
end

function rand_mat(n::Int)
    a = randn(n, n)
    A = Interval.(a)
    mA = MMatrix{n, n}(A)
    sA = SMatrix{n, n}(A)
    return A, mA, sA
end

@testset "Linear Equations" begin

    As = Any[SMatrix{2}(2..3, 1..2, 0..1, 2..3), ]
    bs = Any[SVector(0..120, 60..240), ]
    xs = Any[SVector(-120..90, -60..240), ]

    for i in 2:6
        rand_A = rand_mat(i)[3]
        rand_x = rand_vec(i)[3]
        rand_b = rand_A * rand_x
        push!(As, rand_A)
        push!(bs, rand_b)
        push!(xs, rand_x)
    end

    n = length(As)


    for solver in (gauss_seidel_interval, gauss_seidel_contractor, gauss_elimination_interval)
        @testset "Solver $solver" begin
            #for precondition in (false, true)

                for i in 1:n
                    soln = solver(As[i], bs[i])
                    @test all(xs[i] .⊆ soln)

                end
            #end
        end
    end
end
