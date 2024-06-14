# `roots` interface

## Methods

Three root-finding contractors currently available through the `roots` interface are the following:
  - `Newton` (default);
  - `Krawczyk`;
  - `Bisection`

Both the Newton and Krawczyk methods can determine if a root is unique in an interval, at the cost of requiring that the function is differentiable. The bisection method has no such requirement, but can never guarantee the existence or uniqueness of a root.

The method used is given using the `contractor` keyword argument:

```jl
julia> roots(log, -2..2 ; contractor = Newton)
1-element Vector{Root{Interval{Float64}}}:
 Root([0.999999, 1.00001]_com, :unique)

julia> roots(log, -2..2 ; contractor = Krawczyk)
1-element Vector{Root{Interval{Float64}}}:
 Root([0.999999, 1.00001]_com, :unique)

julia> roots(log, -2..2 ; contractor = Bisection)
1-element Vector{Root{Interval{Float64}}}:
 Root([0.999999, 1.00001]_com, :unknown)
```

Note that as shown in the example, the `log` function does not complain about being given an interval going outside of its domain. While this may be surprising, this is the expected behavior and no root will ever be found outside the domain of a function.

## Explicit derivatives

Newton and Krawczyk methods require the function to be differentiable, but the derivative is usually computed automatically using forward-mode automatic differentiation, provided by the `ForwardDiff.jl` package. It is however possible to provide the derivative explicitly using the `derivative` keyword argument:

```jl
julia> roots(log, -2..2 ; contractor = Newton, derivative = x -> 1/x)
1-element Vector{Root{Interval{Float64}}}:
 Root([0.999999, 1.00001]_com_NG, :unique)

julia> roots(log, -2..2 ; contractor = Krawczyk, derivative = x -> 1/x)
1-element Vector{Root{Interval{Float64}}}:
 Root([0.999999, 1.00001]_com_NG, :unique)
```

When providing the derivative explicitly, the computation is expected to be slightly faster, but the precision of the result is unlikely to be affected.

```jl
julia> using BenchmarkTools

julia> @btime roots(log, -2..2 ; derivative = x -> 1/x)
  4.814 μs (53 allocations: 3.16 KiB)
1-element Vector{Root{Interval{Float64}}}:
 Root([0.999999, 1.00001]_com_NG, :unique)

julia> @btime roots(log, -2..2)
  5.767 μs (53 allocations: 3.16 KiB)
1-element Vector{Root{Interval{Float64}}}:
 Root([0.999999, 1.00001]_com, :unique)
```

This may be useful in some special cases where `ForwardDiff.jl` is unable to compute the derivative of a function. Examples are complex functions and functions whose interval extension must be manually defined (e.g. special functions like `zeta`).

In dimension greater than one, the derivative be given as a function returning the Jacobi matrix:

```jl
julia> function f( (x, y) )
           return [sin(x), cos(y)]
       end
f (generic function with 1 method)

julia> function df( (x, y) )
           return [cos(x) 0 ; 0 -sin(y)]
       end

julia> roots(f, [-3..3, -3..3] ; derivative = df)
2-element Vector{Root{Vector{Interval{Float64}}}}:
 Root(Interval{Float64}[[-1.24409e-21, 1.0588e-22]_com_NG, [-1.5708, -1.57079]_com_NG], :unique)
 Root(Interval{Float64}[[-1.24409e-21, 1.0588e-22]_com_NG, [1.57079, 1.5708]_com_NG], :unique)
```

## Tolerance

An absolute tolerance for the search may be specified as the `abstol` keyword argument.
Currently a method must first be provided in order to be able to choose the tolerance.

```jl
julia> g(x) = sin(exp(x))
g (generic function with 1 method)
        
julia> roots(g, 0..2)
2-element Vector{Root{Interval{Float64}}}:
 Root([1.14472, 1.14474]_com, :unique)
 Root([1.83787, 1.83788]_com, :unique)

julia> roots(g, 0..2 ; abstol = 1e-1)
2-element Vector{Root{Interval{Float64}}}:
 Root([1.14173, 1.15244]_com, :unique)
 Root([1.78757, 1.84273]_com, :unique)
```

A lower tolerance may greatly reduce the computation time, at the cost of an increased number of returned roots having `:unknown` status:

```jl
julia> h(x) = cos(x) * sin(1 / x)
h (generic function with 1 method)

julia> @btime roots(h, 0.05..1)
  79.500 μs (301 allocations: 13.97 KiB)
6-element Vector{Root{Interval{Float64}}}:
 Root([0.0530516, 0.0530517]_com_NG, :unique)
 Root([0.0636619, 0.063662]_com_NG, :unique)
 Root([0.0795774, 0.0795775]_com_NG, :unique)
 Root([0.106103, 0.106104]_com_NG, :unique)
 Root([0.159154, 0.159156]_com_NG, :unique)
 Root([0.318309, 0.31831]_com_NG, :unique)

julia> @btime roots(h, 0.05..1 ; abstol = 1e-2)
  48.500 μs (265 allocations: 11.61 KiB)
6-element Vector{Root{Interval{Float64}}}:
 Root([0.0514445, 0.0531087]_com_NG, :unique)
 Root([0.0570253, 0.0641615]_com, :unknown)
 Root([0.0785458, 0.0797797]_com_NG, :unknown)
 Root([0.104754, 0.107542]_com_NG, :unknown)
 Root([0.157236, 0.165989]_com_NG, :unknown)
 Root([0.31716, 0.319318]_com_NG, :unique)

julia> @btime roots(h, 0.05..1 ; abstol = 1e-1)
  17.300 μs (114 allocations: 5.16 KiB)
3-element Vector{Root{Interval{Float64}}}:
 Root([0.04999999, 0.107542]_com, :unknown)
 Root([0.107541, 0.165989]_com, :unknown)
 Root([0.283803, 0.356099]_com_NG, :unknown)
```

The last example shows a case where the tolerance was too large to be able to isolate the roots in distinct regions.

!!! warning

    For a root `x` of some function, if the absolute tolerance is smaller than `eps(x)` i.e. if `tol + x == x`, `roots` may never be able to converge to the required tolerance and the function may get stuck in an infinite loop.
