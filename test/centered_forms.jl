using IntervalArithmetic, IntervalRootFinding, StaticArrays
using Test



f(x, y) = SVector(x^2 + y^2 - 1, y - 2x)
f(X) = f(X[1], X[2])
X = (-6..6) × (-6..6)

rts = roots(f, X, Bisection)
rts2 = roots(mean_value_form_vector(f), X, Bisection)

@test length(rts) == length(rts2) == 4
