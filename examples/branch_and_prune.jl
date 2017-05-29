
branch_and_prune(-5..5, sin, newton_helper(IntervalRootFinding.K))
branch_and_prune(-5..5, sin, newton_helper(IntervalRootFinding.N))
branch_and_prune(-5..5, sin, bisection_helper)
