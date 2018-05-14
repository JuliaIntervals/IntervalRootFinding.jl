using IntervalArithmetic, StaticArrays, BenchmarkTools, Compat

include("linear_eq.jl")

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
