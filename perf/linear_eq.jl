using IntervalArithmetic, StaticArrays, BenchmarkTools, Compat

include("../src/linear_eq.jl")

function randVec(n::Int)
    a = randn(n)
    A = Interval.(a)
    sA = MVector{n}(A)
    return A, sA
end

function randMat(n::Int)
    a = randn(n, n)
    A = Interval.(a)
    sA = MMatrix{n, n}(A)
    return A, sA
end

function benchmark(max=10)
    for n in 1:max
        A, sA = randMat(n)
        b, sb = randVec(n)
        println("For n = ", n)
        t1 =  @btime gauss_seidel_interval($A, $b)
        println("Array: ", t1)
        t2 =  @btime gauss_seidel_interval_static($sA, $sb)
        println("MArray: ", t2)
        println()
    end
end
