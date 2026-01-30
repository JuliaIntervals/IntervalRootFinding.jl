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

## Tree representation

A branch-and-bound search can be naturally represented as a binary tree: each leaf contains a region and its status and each node represents a bisection. If the tree is known, the topology of the whole search can be reconstructed. There is however a point not determined by the tree.

The branch and bound algorithm used by the `roots` function builds this tree and at the end collect all leaves containing a region with status either `:unknown` or `:unique`. We see below how to access the tree.

## Search strategy

While the tree representation is sufficient to know the topology of the search, it does not determine the order in which leaves are processed during the search. This has an influence on how the tree is created and the amount of memory used, but will not change the resulting tree, unless some limitations on the number of iterations or leaves are enforced.

!!! note

    No such limitations are currently implemented, but they are planned. They will allow to deal, for example, with functions admitting an infinite amount of roots.

Two strategies are currently available: a breadth-first strategy (leaves closer to the root of the tree are processed first); and a depth-first strategy (leaves further away from the root are processed first).

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

## RootProblem and search object

The parameters of a search are represented by a `RootProblem` object
that has the same signature as the `roots` function.

The `RootProblem` can be iterated to perform the search of the roots,
for example to log information at each iteration.

For example, the following stops the search after 15 iterations and
shows the state of the search at that point.

```jldoctest
julia> f(x) = exp(x) - sin(x)
f (generic function with 1 method)

julia> problem = RootProblem(f, interval(-10, 10))
RootProblem{Newton, typeof(f), IntervalRootFinding.var"#9#11"{typeof(f)}, Root{Interval{Float64}}, BranchAndPrune.BreadthFirst, Float64}(Newton, f, IntervalRootFinding.var"#9#11"{typeof(f)}(f), Root([-10.0, 10.0]_com, :unkown), BranchAndPrune.BreadthFirst, 1.0e-7, 0.0, 100000, 0.49609375, true)

julia> state = nothing   # stores current state of the search

julia> for (k, s) in enumerate(problem)
           global state = s
           k == 15 && break  # stop at iteration 15
       end

julia> state.tree
Branching
└─ Branching
   ├─ Branching
   │  ├─ Branching
   │  │  ├─ (:working, Root([-10.0, -8.78861]_com, :unknown)
   │  │  │      Not converged: the root is still being processed)
   │  │  └─ (:working, Root([-8.78861, -7.55814]_com, :unknown)
   │  │         Not converged: the root is still being processed)
   │  └─ (:final, Root([-6.28131, -6.28131]_com, :unique))
   └─ Branching
      ├─ (:working, Root([-5.07782, -3.84735]_com, :unknown)
      │      Not converged: the root is still being processed)
      └─ (:working, Root([-3.84735, -2.5975]_com, :unknown)
             Not converged: the root is still being processed)
```

The elements of the iteration are `SearchState` from the `BranchAndPrune.jl`
package.
In the example, we show the tree that get constructed during the search,
which, at iteration 15, has found one root and have 4 regions of unknown status
to process.