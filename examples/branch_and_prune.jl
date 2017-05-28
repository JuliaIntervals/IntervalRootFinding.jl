branch_and_prune(-1e16..1e16, newton_helper(IntervalRootFinding.K, x->x^2-2)

branch_and_prune(-1e16..1e16, newton_helper(IntervalRootFinding.N, x->x^2-2)

branch_and_prune(-5..5, bisection_helper(sin), 1e-3)
