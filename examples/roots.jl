using IntervalArithmetic
using IntervalArithmetic.Symbols
using IntervalRootFinding

# Find roots of a function in a search region with defaul parameters
rts = roots(sin, -5..5)
rts = roots(sin, -5..6 ; contractor = Bisection, abstol = 1e-1)

# Refine roots with different options
rts = vcat(roots.(sin, rts ; contractor = Krawczyk)...)

# 2D:
f(x, y) = [x^2 + y^2 - 1, y - x]
f(X) = f(X...)  # Function must take a single vector as input

rts = roots(f, [-5..5, -5..5])

# When defining function mixing numbers and intervals use @exact to avoid the NG interval flag
@exact fe(x, y) = [x^2 + y^2 - 1, y - x]
fe(X) = fe(X...)  # Function must take a single vector as input
rts = roots(fe, [-5..5, -5..5])


# From R docs:
# https://www.rdocumentation.org/packages/pracma/versions/1.9.9/topics/broyden
# For performance use SVector
using StaticArrays

function g(x)
    (x1, x2, x3) = x
    SVector(    x1^2 + x2^2 + x3^2 - 1,
                x1^2 + x3^2 - 0.25,
                x1^2 + x2^2 - 4x3
            )
end

X = (-5..5)
rts = roots(g, SVector(X, X, X))

h(xv) = ((x,y) = xv; SVector(2*x - y - exp(-x), -x + 2*y - exp(-y)))
rts = roots(h, SVector(X, X))

# Dennis-Schnabel:
h(xv) = ((x, y) = xv; SVector(x^2 + y^2 - 2, exp(x - 1) + y^3 - 2))
rts = roots(h, SVector(X, X))

# Test suites:
# http://folk.uib.no/ssu029/Pdf_file/Testproblems/testprobRheinboldt03.pdf
# http://www.mat.univie.ac.at/~neum/glopt/test.html
# http://titan.princeton.edu/TestProblems/
# http://www-sop.inria.fr/saga/POL/

##  MINPACK benchmarks: https://github.com/JuliaNLSolvers/NLsolve.jl/blob/master/test/minpack.jl

rosenbrock(xx) = ( (x, y) = xx; SVector( 1 - x, 1000 * (y - x^2) ) )
X = SVector(-1e5..1e5, -1e5..1e5)

rts = roots(rosenbrock, X)

# and other files in NLsolve test suite
# Testing unconstrained optimization software, Moré, Garbow, Hillstrom ACM Trans. Math. Soft. 7 (1), 17-41 (1981)
