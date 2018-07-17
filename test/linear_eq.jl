using IntervalArithmetic, StaticArrays, IntervalRootFinding

function randVec(n::Int)
    a = randn(n)
    A = Interval.(a)
    mA = MVector{n}(A)
    sA = SVector{n}(A)
    return A, mA, sA
end

function randMat(n::Int)
    a = randn(n, n)
    A = Interval.(a)
    mA = MMatrix{n, n}(A)
    sA = SMatrix{n, n}(A)
    return A, mA, sA
end

@testset "Linear Equations" begin

    A = [[2..3 0..1;1..2 2..3], ]
    b = [[0..120, 60..240], ]
    x = [[-120..90, -60..240], ]

    for i in 1:10
        rand_a = randMat(i)[1]
        rand_b = randVec(i)[1]
        rand_c = rand_a * rand_b
        push!(A, rand_a)
        push!(b, rand_c)
        push!(x, rand_b)
    end
    for i in 1:length(A)
        for precondition in (false, true)
            for f in (gauss_seidel_interval, gauss_seidel_contractor, gauss_elimination_interval, \)
                @test all(x[i] .⊆ f(A[i], b[i]))
            end
        end
    end
end
