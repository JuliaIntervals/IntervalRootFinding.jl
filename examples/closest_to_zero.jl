using IntervalArithmetic
using IntervalRootFinding

f(x) = sin(1/(x + 0.5))

pb = RootProblem(f, interval(-10, 10))
state = nothing

for s in pb
    global state = s
    any(isunique, state.roots) && break
end