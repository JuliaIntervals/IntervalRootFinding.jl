using IntervalArithmetic, IntervalRootFinding
using Test

import ForwardDiff

# Rastrigin function:
const A = 10
# n = 2   # dimension of domain

f(x, y) = 2A + x^2 - A*cos(2π*x) + y^2 - A*cos(2π*y)
f(X) = f(X...)

#ForwardDiff.gradient(f, X::IntervalBox) = ForwardDiff.gradient(f, X.v)

# ∇f = X -> ForwardDiff.gradient(f, X)

L = 5.0
X = IntervalBox(-L..(L+1), 2)

rts = IntervalRootFinding.roots(∇(f), X, Newton, 1e-5)
rts2 = IntervalRootFinding.roots(∇(f), X, Krawczyk, 1e-5)

@test length(rts) == length(rts2) == 529
@test all(isunique, rts)
@test all(isunique, rts2)
