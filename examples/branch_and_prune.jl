
branch_and_prune(-5..5, sin, newton_helper(IntervalRootFinding.K))
branch_and_prune(-5..5, sin, newton_helper(IntervalRootFinding.N))
branch_and_prune(-5..5, sin, bisection_helper)

# 2D:
branch_and_prune(   (-5..5) × (-5..5),
                    X -> ( (x,y) = X; IntervalBox(x^2 + y^2 - 1, y - x) ),
                    bisection_helper
                )

branch_and_prune(   (-5..5) × (-5..5),
                    X -> ( (x,y) = X; [x^2 + y^2 - 1, y - x] ),
                    newton_helper(IntervalRootFinding.N)
                )
