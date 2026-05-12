using IntervalArithmetic
using IntervalRootFinding

f(x) = sin(1/(x + 0.5))

## Stop when the first (guaranteed) root is found
pb = RootProblem(f, interval(-10, 10))
state = nothing

for s in pb
    global state = s  # Needed to access the state outside of the loop later
    # `roots(state)` return a list of all roots actually considered
    # Here as soon as one has the :unique status in one of the iteration,
    # we stop the search
    any(isunique, roots(state)) && break
end

# Show the state of the search at the end of the iterations
# Use `roots(state)` to get the list of roots
println(state)

## Choose the order of processing manually
# The `ChangingOrder` search order allows to choose which region
# is processed next at each iteration
pb = RootProblem(f, interval(-10, 10) ; search_order = ChangingOrder)
state = nothing

for s in pb
    global state = s
    any(isunique, roots(state)) && break
    # `unconverged_roots(state)` return the set of roots that have not yet been finalized
    # We need to choose one of them to be processed at the next iteration
    rts = unconverged_roots(state)
    # `mig` return the smallest distance from a member of the interval to zero
    # Therefore the minimal value of `mig` corresponds to the interval closest to zero
    next = argmin(mig.(rts))
    # Use the index to set the next region that should be processed
    set_next!(state.search_order, next)
end

# Show the state of the search at the end of the iterations
# We can see that there are no unprocessed region closest to zero than the found root!
println(state)

# Note that `ChangingOrder` is not optimized,
# defining a new `SearchOrder` may be better in some case