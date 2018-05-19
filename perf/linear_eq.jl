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
    df[:Method] = ["Array", "MArray", "SArray", "Contractor", "ContractorMArray", "ContractorSArray"]
    for n in 1:max
        A, mA, sA = randMat(n)
        b, mb, sb = randVec(n)
        t1 = @belapsed gauss_seidel_interval($A, $b)
        t2 = @belapsed gauss_seidel_interval($mA, $mb)
        t3 = @belapsed gauss_seidel_interval($sA, $sb)
        t4 = @belapsed gauss_seidel_contractor($A, $b)
        t5 = @belapsed gauss_seidel_contractor($mA, $mb)
        t6 = @belapsed gauss_seidel_contractor($sA, $sb)
        df[Symbol("$n")] = [t1, t2, t3, t4, t5, t6]
    end
    a = []
    for i in 1:max
        push!(a, Symbol("$i"))
    end
    df1 = stack(df, a)
    dfnew = unstack(df1, :variable, :Method, :value)
    dfnew = rename(dfnew, :variable => :n)
    println(dfnew)
    dfnew
end

function benchmark_elimination(max=10)
    df = DataFrame()
    df[:Method] = ["Gauss Elimination", "Base.\\"]
    for n in 1:max
        A, mA, sA = randMat(n)
        b, mb, sb = randVec(n)
        t1 = @belapsed gauss_elimination_interval($A, $b)
        t2 = @belapsed gauss_elimination_interval1($A, $b)
        df[Symbol("$n")] = [t1, t2]
    end
    a = []
    for i in 1:max
        push!(a, Symbol("$i"))
    end
    df1 = stack(df, a)
    dfnew = unstack(df1, :variable, :Method, :value)
    dfnew = rename(dfnew, :variable => :n)
    println(dfnew)
    dfnew
end
