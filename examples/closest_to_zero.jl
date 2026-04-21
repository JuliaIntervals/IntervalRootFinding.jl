using IntervalArithmetic
using IntervalRootFinding

f(x) = sin(1/(x + 0.5))

# TODO Explain what is going on here
# Stop when the first (guaranteed) root is found
pb = RootProblem(f, interval(-10, 10))
state = nothing

for s in pb
    global state = s
    any(isunique, roots(state)) && break
end

# Choose the order of processing manually
pb = RootProblem(f, interval(-10, 10) ; search_order = ChangingOrder)
state = nothing

for s in pb
    global state = s
    any(isunique, roots(state)) && break
    state.search_order.next = argmin(mig.(unconverged_roots(state)))
end