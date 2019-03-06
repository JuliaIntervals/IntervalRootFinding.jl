# `roots` interface

## Methods

Three root finding algorithm are currently available through the `roots` interface
  - `Newton` (default)
  - `Krawczyk`
  - `Bisection`

Both Newton and Krawczyk methods can determine if a root is unique in an interval, at the cost of requiring that the function is differentiable. The bisection method as no such requirement, but can never guarantee the existence or unicity of a root.

The method used is given as the third argument of the `roots` function:

```jl
julia> roots(log, -2..2, Newton)
1-element Array{Root{Interval{Float64}},1}:
 Root([0.999996, 1.00001], :unique)

julia> roots(log, -2..2, Krawczyk)
1-element Array{Root{Interval{Float64}},1}:
 Root([0.999984, 1.00002], :unique)

julia> roots(log, -2..2, Bisection)
1-element Array{Root{Interval{Float64}},1}:
 Root([0.999454, 1.00039], :unknown)
```

Note that as shown in the example, the `log` function does not complain about being given an interval going outside of its domain. While maybe surprising, this is the expected behavior and no root will ever be found outside the domain of a function.

## Explicit derivatives

Newton and Krawczyk methods require the function to be differentiable but the derivative is computed automatically using `ForwardDiff.jl`. It is however possible to provide the derivative explicitly for these methods as the second argument of the `roots` function:

```jl
julia> roots(log, x -> 1/x, -2..2, Newton)
1-element Array{Root{Interval{Float64}},1}:
 Root([0.999996, 1.00001], :unique)

julia> roots(log, x -> 1/x, -2..2, Krawczyk)
1-element Array{Root{Interval{Float64}},1}:
 Root([0.999984, 1.00002], :unique)
```

When providing the derivative explicitly, the computation is expected to be slightly faster, but the precision of the result is unlikely to be affected.

```jl
julia> using BenchmarkTools

julia> @btime roots(log, x -> 1/x, -2..2, Newton)
  38.600 μs (371 allocations: 27.01 KiB)
1-element Array{Root{Interval{Float64}},1}:
 Root([0.999996, 1.00001], :unique)

julia> @btime roots(log, -2..2, Newton)
  51.799 μs (373 allocations: 27.20 KiB)
1-element Array{Root{Interval{Float64}},1}:
 Root([0.999996, 1.00001], :unique)
```

It may be useful in some special case where `ForwardDiff.jl` fail to compute the derivative of a function. Currently the only known case of this happening is for functions whose interval extension must be manually define (e.g. special functions like `zeta`).

In dimension greater than one, the function of interest must return a `SVector`, a type provided by the `StaticArrays` package, but otherwise work in the same way as in the 1D case.

```jl
julia> function f(xx)
           x, y = xx
           return SVector(sin(x), cos(y))
       end
f (generic function with 1 method)

julia> roots(f, IntervalBox(-3..3, 2))
2-element Array{Root{IntervalBox{2,Float64}},1}:
 Root([-1.13556e-19, 2.3664e-20] × [1.57079, 1.5708], :unique)
 Root([-7.92188e-19, 3.20973e-19] × [-1.5708, -1.57079], :unique)
```

The corresponding to the 1D derivative is a function returning the Jacobi matrix of the function as a `SMatrix`:

```jl
julia> function df(xx)
           x, y = xx
           return SMatrix{2, 2}(cos(x), 0, 0, -sin(y))
       end
df (generic function with 1 method)

julia> roots(f, df, IntervalBox(-3..3, 2), Newton)
2-element Array{Root{IntervalBox{2,Float64}},1}:
 Root([-2.35877e-07, 8.22858e-07] × [1.57079, 1.5708], :unique)
 Root([-7.19393e-07, 1.55473e-06] × [-1.5708, -1.57079], :unique)
```


## Tolerance

Absolute tolerance for the search can be specified as the last argument of the `roots` function, default being `1e-15`. Currently, a method must be provided to be able to choose the tolerance.

```jl
julia> g(x) = sin(exp(x))
g (generic function with 1 method)

julia> roots(g, 0..2, Newton)
2-element Array{Root{Interval{Float64}},1}:
 Root([1.83787, 1.83788], :unique)
 Root([1.14472, 1.14473], :unique)

julia> roots(g, 0..2, Newton, 1e-2)
2-element Array{Root{Interval{Float64}},1}:
 Root([1.83745, 1.83974], :unique)
 Root([1.14471, 1.14475], :unique)
```

Lower tolerance may greatly reduce computation time at the cost of increased number of returned roots having `:unkown` status.

```jl
julia> h(x) = cos(x)*sin(1/x)
h (generic function with 1 method)

julia> @btime roots(h, 0.05..1, Newton)
  1.316 ms (9676 allocations: 202.21 KiB)
6-element Array{Root{Interval{Float64}},1}:
 Root([0.106103, 0.106104], :unique)
 Root([0.318309, 0.31831], :unique)
 Root([0.0795774, 0.0795775], :unique)
 Root([0.0636619, 0.063662], :unique)
 Root([0.0530516, 0.0530517], :unique)
 Root([0.159154, 0.159155], :unique)             

julia> @btime roots(h, 0.05..1, Newton, 1e-2)
   475.500 μs (4171 allocations: 94.00 KiB)
 6-element Array{Root{Interval{Float64}},1}:
  Root([0.317179, 0.319299], :unique)
  Root([0.157209, 0.165989], :unknown)
  Root([0.104739, 0.107542], :unknown)
  Root([0.0515382, 0.0531049], :unique)
  Root([0.0785458, 0.0797755], :unknown)
  Root([0.0570253, 0.0641614], :unknown)

julia> @btime roots(h, 0.05..1, Newton, 1e-1)
  151.101 μs (1382 allocations: 32.86 KiB)
3-element Array{Root{Interval{Float64}},1}:
 Root([0.107541, 0.165989], :unknown)
 Root([0.283803, 0.3555], :unknown)
 Root([0.0499999, 0.107542], :unknown)     
```

The last example show a tolerance too big to be able to isolate the roots in distinct regions.

!!! warning

    For a root `x` of some function, if the absolute tolerance is smaller than `eps(x)` i.e. if `tol + x == x`, `roots` may never be able to converge to the require tolerance and the function may get stuck in an infinite loop.
