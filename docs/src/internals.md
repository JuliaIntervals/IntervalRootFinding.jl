```@meta
DocTestSetup = quote
    using IntervalArithmetic, IntervalArithmetic.Symbols, IntervalRootFinding
end
```
# Internals

This section describes some of the internal mechanism of the package and several ways to use them to customize a search.

## Branch and bound

When `roots` is called, it performs a **branch-and-bound** search, that is, it iteratively looks at a region $X$ and for each region tries to determine if it contains a root. It then does the following:
  - If $X$ is proven to contain *no* root, it discards it.
  - If $X$ is proven to contain *exactly one root*, it tries to get the best possible bounds for the root and store the resulting region in the list of `Root`s to output with a `:unique` status.
  - If the test is inconclusive and the size of $X$ is smaller than the tolerance, it stores it in the list of `Root`s to output with `:unknown` status.
  - If the test is inconclusive and the size of $X$ is larger than the tolerance, it bisects $X$ and then processes each resulting half.

At some point all regions will either have a determined status or be smaller than the tolerance, and the algorithm will halt and return all stored roots.

## Contractors

To determine the status of a region, the algorithm uses so-called *contractors*.
When we contract a region (wrapped in a `Root` object),
it returns the status of the root and the region.
The contractors are various methods to guarantee and refine the
status of a root.
The available contractors are `Bisection`, `Newton` or `Krawczyk`.

```jldoctest
julia> using IntervalRootFinding: contract

julia> contract(Newton, sin, cos, Root(pi ± 0.001, :unknown))
Root([3.14159, 3.14159]_com, :unique)

julia> contract(Newton, sin, cos, Root(2 ± 0.001, :unknown))
Root([1.999, 2.001]_com, :empty)
```

While it is the fastest and doesn't require derivatives,
`Bisection` can never prove the existence or unicity of a root.

## Tree representation

A branch-and-bound search can be naturally represented as a binary tree: each leaf contains a region and its status and each node represents a bisection. If the tree is known, the topology of the whole search can be reconstructed. There is however a point not determined by the tree.

The branch and bound algorithm used by the `roots` function builds this tree
and at the end collect all leaves containing a region with status either `:unknown` or `:unique`.

## Search strategy

While the tree representation is sufficient to know the topology of the search,
it does not determine the order in which leaves are processed during the search.
This has an influence on how the tree is created and the amount of memory used,
but will not change the resulting tree,
unless the maximal number of iterations is reached.

Two strategies are currently available: a breadth-first strategy (process all regions before processing a sub-region);
and a depth-first strategy (immediately process the sub-regions of the last processed region).

Usually, the breadth-first strategy (default) is preferred,
as the depth-first strategy may get "stuck" on a single point with a singularity,
until either the tolerance or the maximum number of iterations is reached.

## RootProblem and search object

The parameters of a search are represented by a `RootProblem` object
that has the same signature as the `roots` function.

The `RootProblem` can be iterated to perform the search of the roots,
to log information, introduce an early stop criterion,
or even modify the iteration state.

The follow example stops the search after 15 iterations and
shows the state of the search at that point.

```jldoctest
julia> f(x) = exp(x) - sin(x)
f (generic function with 1 method)

julia> problem = RootProblem(f, interval(-10, 10))
RootProblem
  Contractor: Newton
  Function: f
  Search region: [-10.0, 10.0]_com
  Search order: BreadthFirst
  Absolute tolerance: 1.0e-7
  Relative tolerance: 0.0
  Maximum iterations: 100000
  Ignored errors: DataType[IntervalArithmetic.InconclusiveBooleanOperation]

julia> state = nothing   # stores current state of the search, so that it will be available outside of the loop

julia> for s in problem
           global state = s
           state.iteration == 15 && break
       end

julia> state
SearchState{BreadthFirst{Root{Interval{Float64}}}, Root{Interval{Float64}}}
  iteration: 15
  5 regions being processed
    Root([-0.078125, 1.15234]_com, :unknown)
    Root([-3.84735, -2.5975]_com, :unknown)
    Root([-5.07782, -3.84735]_com, :unknown)
    Root([-8.78861, -7.55814]_com, :unknown)
    Root([-10.0, -8.78861]_com, :unknown)
  1 finalized regions
    Root([-6.28131, -6.28131]_com, :unique)

```

The elements of the iteration are `SearchState` from the `BranchAndPrune.jl` package.

The roots being processed can be accessed with `unfinished_roots(state)`
and the finalized one with `finished_roots(state)`,
while `roots(state)` return a list containing both.

See `examples/closest_to_zero.jl` for a more advanced example,
which shows how you can find the root that is closest to zero
and immediately terminate the search.