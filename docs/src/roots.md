```@meta
DocTestSetup = quote
    using IntervalArithmetic, IntervalArithmetic.Symbols, IntervalRootFinding
end
```
# `roots` interface

## Methods

Three root-finding contractors currently available through the `roots` interface are the following:
  - `Newton` (default);
  - `Krawczyk`;
  - `Bisection`

Both the Newton and Krawczyk methods can determine if a root is unique in an interval, at the cost of requiring that the function is differentiable. The bisection method has no such requirement, but can never guarantee the existence or uniqueness of a root.

The method used is given using the `contractor` keyword argument:

```jldoctest
julia> roots(log, -2..2 ; contractor = Newton)
1-element Vector{Root{Interval{Float64}}}:
 Root([1.0, 1.0]_com, :unique)

julia> roots(log, -2..2 ; contractor = Krawczyk)
1-element Vector{Root{Interval{Float64}}}:
 Root([1.0, 1.0]_com, :unique)

julia> roots(log, -2..2 ; contractor = Bisection)
1-element Vector{Root{Interval{Float64}}}:
 Root([1.0, 1.0]_com, :unknown)
    Not converged: region size smaller than the tolerance
```

Note that as shown in the example, the `log` function does not complain about being given an interval going outside of its domain. While this may be surprising, this is the expected behavior and no root will ever be found outside the domain of a function.

## Explicit derivatives

Newton and Krawczyk methods require the function to be differentiable, but the derivative is usually computed automatically using forward-mode automatic differentiation, provided by the `ForwardDiff.jl` package. It is however possible to provide the derivative explicitly using the `derivative` keyword argument:

```jldoctest
julia> roots(log, -2..2 ; contractor = Newton, derivative = x -> 1/x)
1-element Vector{Root{Interval{Float64}}}:
 Root([1.0, 1.0]_com_NG, :unique)

julia> roots(log, -2..2 ; contractor = Krawczyk, derivative = x -> 1/x)
1-element Vector{Root{Interval{Float64}}}:
 Root([1.0, 1.0]_com_NG, :unique)
```

When providing the derivative explicitly, the computation is expected to be slightly faster, but the precision of the result is unlikely to be affected.

```julia-repl
julia> using BenchmarkTools

julia> @btime roots(log, -2..2 ; derivative = x -> 1/x)
  7.050 μs (129 allocations: 9.33 KiB)
1-element Vector{Root{Interval{Float64}}}:
 Root([1.0, 1.0]_com_NG, :unique)

julia> @btime roots(log, -2..2)
  7.743 μs (129 allocations: 9.33 KiB)
1-element Vector{Root{Interval{Float64}}}:
 Root([1.0, 1.0]_com, :unique)
```

This may be useful in some special cases where `ForwardDiff.jl` is unable to compute the derivative of a function. Examples are complex functions and functions whose interval extension must be manually defined (e.g. special functions like `zeta`).

In dimension greater than one, the derivative can be given as a function returning the Jacobi matrix:

```jldoctest
julia> f( (x, y) ) = [sin(x), cos(y)]
f (generic function with 1 method)

julia> df( (x, y) ) = [cos(x) 0 ; 0 -sin(y)]
df (generic function with 1 method)

julia> roots(f, [-3..3, -3..3] ; derivative = df)
2-element Vector{Root{Vector{Interval{Float64}}}}:
 Root(Interval{Float64}[[-1.28378e-21, 1.58819e-22]_com_NG, [-1.5708, -1.5708]_com_NG], :unique)
 Root(Interval{Float64}[[-1.28378e-21, 1.58819e-22]_com_NG, [1.5708, 1.5708]_com_NG], :unique)
```

## Tolerance

An absolute tolerance for the search may be specified as the `abstol` keyword argument.

```jldoctest
julia> g(x) = sin(exp(x))
g (generic function with 1 method)

julia> roots(g, 0..2)
2-element Vector{Root{Interval{Float64}}}:
 Root([1.14473, 1.14473]_com, :unique)
 Root([1.83788, 1.83788]_com, :unique)

julia> roots(g, 0..2 ; abstol = 1e-1)
2-element Vector{Root{Interval{Float64}}}:
 Root([1.14173, 1.15244]_com, :unique)
 Root([1.78757, 1.84272]_com, :unique)
```

A lower tolerance may greatly reduce the computation time, at the cost of an increased number of returned roots having `:unknown` status:

```julia-repl
julia> h(x) = cos(x) * sin(1 / x)
h (generic function with 1 method)

julia> @btime roots(h, 0.05..1)
  484.488 μs (12218 allocations: 590.76 KiB)
6-element Vector{Root{Interval{Float64}}}:
 Root([0.0530516, 0.0530517]_com_NG, :unique)
 Root([0.063662, 0.063662]_com_NG, :unique)
 Root([0.0795775, 0.0795775]_com_NG, :unique)
 Root([0.106103, 0.106103]_com_NG, :unique)
 Root([0.159155, 0.159155]_com_NG, :unique)
 Root([0.31831, 0.31831]_com_NG, :unique)

julia> @btime roots(h, 0.05..1 ; abstol = 1e-2)
  249.720 μs (6141 allocations: 297.80 KiB)
6-element Vector{Root{Interval{Float64}}}:
 Root([0.0514446, 0.0531086]_com_NG, :unique)
 Root([0.0570254, 0.0641614]_com, :unknown)
    Not converged: region size smaller than the tolerance
 Root([0.0785458, 0.0797797]_com_NG, :unknown)
    Not converged: region size smaller than the tolerance
 Root([0.104755, 0.107541]_com_NG, :unknown)
    Not converged: region size smaller than the tolerance
 Root([0.157236, 0.165988]_com_NG, :unknown)
    Not converged: region size smaller than the tolerance
 Root([0.317161, 0.319318]_com_NG, :unique)

julia> @btime roots(h, 0.05..1 ; abstol = 1e-1)
  66.330 μs (1402 allocations: 69.10 KiB)
3-element Vector{Root{Interval{Float64}}}:
 Root([0.05, 0.107541]_com, :unknown)
    Not converged: region size smaller than the tolerance
 Root([0.107541, 0.165988]_com, :unknown)
    Not converged: region size smaller than the tolerance
 Root([0.283804, 0.356099]_com_NG, :unknown)
    Not converged: region size smaller than the tolerance
```

The last example shows a case where the tolerance was too large to be able to isolate the roots in distinct regions.

!!! warning

    For a root `x` of some function, if the absolute tolerance is smaller than `eps(x)` i.e. if `tol + x == x`, `roots` may never be able to converge to the required tolerance and the function may get stuck in an infinite loop.
    To avoid that, either increase `abstol` or `reltol`,
    or set a maximum number of iterations with the `max_iteration` keyword.

## Error during computation

By default, the `roots` function will ignore errors
and continue bisecting the given interval,
hoping to eventually find regions that are small enough to avoid the error.

This is particularly useful when using code that is using comparisons,
which are ill-defined for intervals
(see [the IntervalArithmetic.jl documentation](https://juliaintervals.github.io/IntervalArithmetic.jl/stable/manual/usage/#Comparisons)
for more information).

For example, the following will find all roots,
and flag the erroring singularity:
```jldoctest bisect_on_error
julia> f(x) = x > 2 ? 7 - x : x
f (generic function with 1 method)

julia> roots(f, -10 .. 10)
3-element Vector{Root{Interval{Float64}}}:
 Root([0.0, 0.0]_com, :unique)
 Root([2.0, 2.0]_com, :unknown)
    Not converged: region size smaller than the tolerance
    Warning: an error was encountered during computation
 Root([7.0, 7.0]_com_NG, :unique)
```

This behavior can be deactivated by setting `bisect_on_error` to false,
interrupting the process as soon as an error is encountered.
This is notably useful while debugging your code.

```julia-repl
julia> roots(f, -10 .. 10 ; bisect_on_error = false)
ERROR: ArgumentError: `<` is purposely not supported for overlapping intervals. See instead `strictprecedes`
```

!!! warning

    This behavior is considered experimental and subject to change.