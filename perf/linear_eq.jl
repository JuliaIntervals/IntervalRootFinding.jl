using IntervalArithmetic, StaticArrays, BenchmarkTools, Compat

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

function benchmark(max=10)
    for n in 1:max
        A, mA, sA = randMat(n)
        b, mb, sb = randVec(n)
        println("For n = ", n)
        t1 =  @btime gauss_seidel_interval($A, $b)
        println("Array: ", t1)
        t2 =  @btime gauss_seidel_interval_static($mA, $mb)
        println("MArray: ", t2)
        t3 =  @btime gauss_seidel_interval_static1($sA, $sb)
        println("SArray1: ", t3)
        t4 =  @btime gauss_seidel_interval_static2($sA, $sb)
        println("SArray2: ", t4)
        println()
    end
end
