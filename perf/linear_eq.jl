using BenchmarkTools, Compat, DataFrames, IntervalRootFinding, IntervalArithmetic, StaticArrays

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
    df = DataFrame()
    df[:Method] = ["Array", "MArray", "SArray1", "SArray2", "Contractor"]
    for n in 1:max
        A, mA, sA = randMat(n)
        b, mb, sb = randVec(n)
        t1 = @belapsed gauss_seidel_interval($A, $b)
        t2 = @belapsed gauss_seidel_interval_static($mA, $mb)
        t3 = @belapsed gauss_seidel_interval_static1($sA, $sb)
        t4 = @belapsed gauss_seidel_interval_static2($sA, $sb)
        t5 = @belapsed gauss_seidel_contractor($A, $b)
        df[Symbol("n = $n")] = [t1, t2, t3, t4, t5]
    end
    println(df)
    df
end
