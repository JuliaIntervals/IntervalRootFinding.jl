
using IntervalArithmetic, IntervalRootFinding, StaticArrays

roots = branch_and_prune(-5..5, sin, Bisection)

roots = branch_and_prune(roots, sin, Newton)

roots = branch_and_prune(-5..6, sin, Bisection)
roots = branch_and_prune(roots, sin, Newton)


# 2D:
roots = branch_and_prune(   (-5..5) × (-5..5),
                            X -> ( (x,y) = X; IntervalBox(x^2 + y^2 - 1, y - x) ),
                            Bisection
                        )

roots = branch_and_prune(   (-5..5) × (-5..5),
                            X -> ( (x,y) = X; SVector(x^2 + y^2 - 1, y - x) )
                            Newton
                        )

roots = branch_and_prune(   roots,
                            X -> ( (x,y) = X; SVector(x^2 + y^2 - 1, y - x) ),
                            Newton
                        )

# complex:
x = -5..6
Xc = x + im * x
f(z) = z^3 - 1

roots = branch_and_prune(Xc, f, Bisection)
roots = branch_and_prune(roots, f, Newton)


# From R docs:

# https://www.rdocumentation.org/packages/pracma/versions/1.9.9/topics/broyden

function f(x)
    x1, x2, x3 = x
    SVector(    x1^2 + x2^2 + x3^2 - 1,
                x1^2 + x3^2 - 0.25,
                x1^2 + x2^2 - 4x3
            )
end

X = (-5..5)
@time roots = branch_and_prune(X × X × X, f, Bisection)
@time roots = branch_and_prune(roots, f, Newton)
length(roots)





g(x) = SVector(2*x[1] - x[2] - exp(-x[1]), -x[1] + 2*x[2] - exp(-x[2]))

roots = branch_and_prune(X × X, g, Bisection)
roots = branch_and_prune(roots, g, Newton)



# Dennis-Schnabel:
h(x) = SVector(x[1]^2 + x[2]^2 - 2, exp(x[1] - 1) + x[2]^3 - 2)

roots = branch_and_prune(X × X, h, Bisection)
roots = branch_and_prune(roots, h, Newton)

# Test suites:

# http://folk.uib.no/ssu029/Pdf_file/Testproblems/testprobRheinboldt03.pdf

# http://www.mat.univie.ac.at/~neum/glopt/test.html

# http://titan.princeton.edu/TestProblems/

# http://www-sop.inria.fr/saga/POL/

# The MINPACK benchmarks () https://github.com/JuliaNLSolvers/NLsolve.jl/blob/master/test/minpack.jl
# and other files in NLsolve test suite

# Testing unconstrained optimization software, Moré, Garbow, Hillstrom ACM Trans. Math. Soft. 7 (1), 17-41 (1981)
