
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
