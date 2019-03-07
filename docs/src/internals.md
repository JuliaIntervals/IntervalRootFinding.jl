# Internals

This section describe some of the internal mechanism of the package and several ways to use them to customize a search.

## Branch and bound

When `roots` is called, it performed a branch and bound search, that is it iteratively look at region and for each region try to determine if it contains a root. It then do the following
  - If the region is proven to contain *no* root, it discards it.
  - If the region is proven to contain *exactly one root*, it tries to get the best possible bounds for the root and store the resulting region in the list of `Root`s to output with `:unique` status.
  - If the test is inconclusive and the size of the region smaller than the tolerance, store it in the list of `Root`s to output with `:unkown` status.
  - If the test is inconclusive and the size of the region bigger than the tolerance, it bisects the region and then process each of resulting half.

At some point all regions will either have a determined status or be smaller than the tolerance, and the algorithm will halt and return all stored roots.

## Tree representation

A branch an bound search can be naturally represented as a binary tree: each leaf contains a region and its status and each node represents a bisection. If the tree is known, the topology of the whole search can be reconstructed. There is however a point not determined by the tree.

The branch and bound algorithm used by the `roots` function builds this tree and at the end collect all leaves containing a region with status either `:unknown` or `:unique`. We see below how to access the tree.

## Search strategy

While the tree representation is sufficient to know the topology of the search, it does not determine the order in which leaves are process during the search. This has an influence on how the tree is created and the amount of memory used, but will not change the resulting tree unless some limitations on the number of iterations or leaves are enforced.

!!! note

    No such limitations are currently implemented but they are planned. They will allow to deal, for example, with functions admitting an infinite amount of roots.

Two strategy are currently available, a breadth first strategy (the leaves closer to the root of the tree are processed first) and a depth first strategy (the leaves further away from the root are processed first).

## Contractors

To determine the status of a region, the algorithm use so called *contractors*. Contractor a callable object build from a function (`Bisection`) and sometimes its derivative (`Newton` and `Krawczyk`). When called with on a region (wrapped in a `Root` object) and a tolerance it returns the status of the root and the region (refined if the region contained an unique root).

```jl
julia> C = Newton(sin, cos)
Newton{typeof(sin),typeof(cos)}(sin, cos)

julia> C(Root(pi ± 0.001, :unknown), 1e-10)
Root([3.14159, 3.1416], :unique)

julia> C(Root(2 ± 0.001, :unkown), 1e-10)
Root([1.99899, 2.00101], :empty)
```

Contractors play a central role in the algorithm. They are indeed the only part of it that varies for different methods.

## Search object

Now that we have presented the foundation of the internal algorithm, we can discuss the representation of the search. Each search strategy has a type associated, the defined types being `BreadthFirstSearch` and `DepthFirstSearch`.

A search must be given three information:
  1. The region to search.
  2. A contractor.
  3. A tolerance.

```jl
julia> f(x) = exp(x) - sin(x)
f (generic function with 1 method)

julia> df(x) = exp(x) - cos(x)
df (generic function with 1 method)

julia> C = Newton(f, df)
Newton{typeof(f),typeof(df)}(f, df)

julia> search = DepthFirstSearch(-10..10, C, 1e-10)
DepthFirstSearch{Interval{Float64},Newton{typeof(f),typeof(df)},
Float64}(Root([-10, 10], :unknown), Newton{typeof(f),typeof(df)}(f, df),
1.0e-10)
```

Then the search is performed using the iterator interface, i.e. a for loop.

```jl
julia> endtree = nothing

julia> for tree in search
           global endtree = tree
       end

julia> endtree
Working tree with 9 elements of type Root{Interval{Float64}}
Indices: [1, 2, 3, 4, 5, 7, 8, 9, 10]
Structure:
  [1] Node with children [2]
    [2] Node with children [3, 4]
      [3] Node with children [8, 9]
        [8] Node with children [10]
          [10] Leaf (:final) with data Root([-9.42486, -9.42485], :unique)
        [9] Leaf (:final) with data Root([-6.28132, -6.28131], :unique)
      [4] Node with children [5]
        [5] Node with children [7]
          [7] Leaf (:final) with data Root([-3.18307, -3.18306], :unique)
```

The element of the iteration are the trees (of type `BBTree`) that get constructed during the search. In the above example we simply get the final iteration of the tree and show it. The list of final roots can be retrieved using the `data` function:

```jl
julia> data(endtree)
3-element Array{Root{Interval{Float64}},1}:
 Root([-3.18307, -3.18306], :unique)
 Root([-6.28132, -6.28131], :unique)
 Root([-9.42486, -9.42485], :unique)
```

We can use this interface to do some analysis at each iteration.

```jl
julia> for (k, tree) in enumerate(search)
           println("The tree at iteration $k has $(IntervalRootFinding.nnodes(tree)) nodes")
       end
The tree at iteration 1 has 3 nodes
The tree at iteration 2 has 5 nodes
The tree at iteration 3 has 4 nodes
                ⋮   # lines not shown to keep it short
The tree at iteration 17 has 10 nodes
The tree at iteration 18 has 9 nodes
The tree at iteration 19 has 9 nodes
```
