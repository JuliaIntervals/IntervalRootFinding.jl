
using IntervalArithmetic, IntervalRootFinding, StaticArrays

rts = roots(sin, -5..5)

rts = roots(sin, -5..6, Bisection)
rts = roots(sin, rts, Newton)

# 2D:
f(x, y) = SVector(x^2 + y^2 - 1, y - x)
f(X) = f(X...)

rts = roots(f, (-5..5) × (-5..5))
rts = roots(f, rts, Bisection)

# complex:
x = -5..6
Xc = x + im * x
f(z) = z^3 - 1

rts = roots(Xc, f, Bisection)
rts = roots(rts, f, Newton)


# From R docs:

# https://www.rdocumentation.org/packages/pracma/versions/1.9.9/topics/broyden

function g(x)
    x1, x2, x3 = x
    SVector(    x1^2 + x2^2 + x3^2 - 1,
                x1^2 + x3^2 - 0.25,
                x1^2 + x2^2 - 4x3
            )
end

X = (-5..5)
@time rts = roots(g, X × X × X)





h(x) = SVector(2*x[1] - x[2] - exp(-x[1]), -x[1] + 2*x[2] - exp(-x[2]))

rts = roots(h, X × X, Bisection)
rts = roots(rts, g, Newton)



# Dennis-Schnabel:
h(x) = SVector(x[1]^2 + x[2]^2 - 2, exp(x[1] - 1) + x[2]^3 - 2)

rts = roots(X × X, h, Bisection)
rts = roots(rts, h, Newton)

# Test suites:

# http://folk.uib.no/ssu029/Pdf_file/Testproblems/testprobRheinboldt03.pdf

# http://www.mat.univie.ac.at/~neum/glopt/test.html

# http://titan.princeton.edu/TestProblems/

# http://www-sop.inria.fr/saga/POL/

# The MINPACK benchmarks () https://github.com/JuliaNLSolvers/NLsolve.jl/blob/master/test/minpack.jl
# and other files in NLsolve test suite

# Testing unconstrained optimization software, Moré, Garbow, Hillstrom ACM Trans. Math. Soft. 7 (1), 17-41 (1981)
