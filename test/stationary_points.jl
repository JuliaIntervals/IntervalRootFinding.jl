using IntervalArithmetic, IntervalRootFinding
using Base.Test

import ForwardDiff

# Rastrigin function:
const A = 10
# n = 2   # dimension of domain

f(x, y) = 2A + x^2 - A*cos(2π*x) + y^2 - A*cos(2π*y)
f(X) = f(X...)

ForwardDiff.gradient(f, X::IntervalBox) = ForwardDiff.gradient(f, X.v)

∇f = X -> ForwardDiff.gradient(f, X)

L = 5.0
X = IntervalBox( (-L..L+1), 2)

@time rts = IntervalRootFinding.roots(∇f, X, Bisection)
@time rts = IntervalRootFinding.roots(∇f, rts, Newton, 1e-20)

@test length(rts) == 529
