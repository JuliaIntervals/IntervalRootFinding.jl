# IntervalRootFinding.jl

[![Build Status](https://travis-ci.org/JuliaIntervals/IntervalRootFinding.jl.svg?branch=master)](https://travis-ci.org/JuliaIntervals/IntervalRootFinding.jl)

[![codecov.io](http://codecov.io/github/JuliaIntervals/IntervalRootFinding/coverage.svg?branch=master)](http://codecov.io/github/JuliaIntervals/IntervalRootFinding.jl?branch=master)

This package provides guaranteed methods for finding **roots** of functions $f: \mathbb{R}^n \to \mathbb{R}^n$ with $n \ge 1$, i.e. vectors (or scalars, for $n=1$) $\mathbb{x}$ for which $f(\mathbb{x}) = \mathbb{0}$. It guarantees to find *all* roots inside a given box in $\mathbb{R}^n$, or report subboxes for which it is unable to provide guarantees.

To do so, it uses methods from interval analysis, using interval arithmetic from the [`IntervalArithmetic.jl`](https://github.com/JuliaIntervals/IntervalArithmetic.jl) package by the same authors.


## Basic usage examples

The basic function is `roots`. A standard Julia function and a box (interval in 1D) are supplied. Optional search methods (currently `Bisection` or `Newton`) and tolerances may be provided.

### 1D

```jl
julia> using IntervalArithmetic, IntervalRootFinding

julia> rts = roots(x->x^2 - 2, -10..10, Bisection)
2-element Array{IntervalRootFinding.Root{IntervalArithmetic.Interval{Float64}},1}:
 Root([1.41418, 1.4148], :unknown)
 Root([-1.4148, -1.41418], :unknown)
```
An array is returned that consists of `Root` objects, containing an interval and the status of that interval. Here, `:unknown` indicates that there *may* be a root in the interval. Any region not contained in one of the intervals is guaranteed *not* to contain a root.

```jl
julia> rts = roots(x->x^2 - 2, -10..10)   # default is Newton
2-element Array{IntervalRootFinding.Root{IntervalArithmetic.Interval{Float64}},1}:
 Root([1.41421, 1.41422], :unique)
 Root([-1.41422, -1.41421], :unique)
```
The interval Newton method used here can *guarantee* that there exists a unique root in each of these intervals. Again, other regions have been excluded.

Interval methods are not able to control multiple roots:
```jl
julia> g(x) = (x^2-2)^2 * (x^2 - 3)
g (generic function with 1 method)

julia> roots(g, -10..10)
4-element Array{IntervalRootFinding.Root{IntervalArithmetic.Interval{Float64}},1}:
 Root([1.73205, 1.73206], :unique)
 Root([1.41418, 1.4148], :unknown)
 Root([-1.4148, -1.41418], :unknown)
 Root([-1.73206, -1.73205], :unique)
 ```
 The two double roots are reported as being possible roots, but no guarantees are given. The single roots are guaranteed to exist and be unique within the corresponding intervals.
 ```


### nD

For dimensions $n>1$, functions must return a Julia vector or an `SVector` from the `StaticArrays.jl` package.

The Rosenbrock function is well known to be difficult to optimize, but is not a problem for this package:

```jl
julia> using StaticArrays;

julia> rosenbrock(xx) = ( (x, y) = xx; SVector( 1 - x, 100 * (y - x^2) ) );

julia> X = IntervalBox(-1e5..1e5, 2)  # 2D IntervalBox;

julia> rts = roots(rosenbrock, X)
1-element Array{IntervalRootFinding.Root{IntervalArithmetic.IntervalBox{2,Float64}},1}:
 Root([1, 1] × [1, 1], :unique)
 ```
 Again, a unique root has been found.


A 3D example:
```jl
julia> function g(x)
           (x1, x2, x3) = x
           SVector(    x1^2 + x2^2 + x3^2 - 1,
                       x1^2 + x3^2 - 0.25,
                       x1^2 + x2^2 - 4x3
                   )
       end
g (generic function with 1 method)

julia> X = (-5..5)
[-5, 5]

julia> @time rts = roots(g, X × X × X)
  0.843263 seconds (766.37 k allocations: 36.338 MiB, 1.12% gc time)
4-element Array{IntervalRootFinding.Root{IntervalArithmetic.IntervalBox{3,Float64}},1}:
 Root([0.440762, 0.440763] × [0.866025, 0.866026] × [0.236067, 0.236068], :unique)
 Root([0.440762, 0.440763] × [-0.866026, -0.866025] × [0.236067, 0.236068], :unique)
 Root([-0.440763, -0.440762] × [0.866025, 0.866026] × [0.236067, 0.236068], :unique)
 Root([-0.440763, -0.440762] × [-0.866026, -0.866025] × [0.236067, 0.236068], :unique)
 ```
 There are guaranteed to be four unique roots.

### Stationary points

Stationary points of a function $f:\mathbb{R}^n \to \mathbb{R}$ may be found as zeros of the gradient.
The package exports the `∇` operator to calculate gradients using `ForwardDiff.jl`:

```jl
julia> f(xx) = ( (x, y) = xx; sin(x) * sin(y) )
f (generic function with 1 method)

julia> ∇f = ∇(f)  # gradient operator from the package
(::#53) (generic function with 1 method)

julia> rts = roots(∇f, IntervalBox(-5..6, 2), Newton, 1e-5)
25-element Array{IntervalRootFinding.Root{IntervalArithmetic.IntervalBox{2,Float64}},1}:
 Root([4.71238, 4.71239] × [4.71238, 4.71239], :unique)
 Root([4.71238, 4.71239] × [1.57079, 1.5708], :unique)
 ⋮
 [output snipped for brevity]
```



Some basic documentation (mostly now out of date) is available at [here](https://juliaintervals.github.io/IntervalRootFinding.jl/latest/).

## Authors
- [Luis Benet](http://www.cicc.unam.mx/~benet/), Instituto de Ciencias Físicas,
Universidad Nacional Autónoma de México (UNAM)
- [David P. Sanders](http://sistemas.fciencias.unam.mx/~dsanders),
Departamento de Física, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM)

## Acknowledgements ##

Financial support is acknowledged from DGAPA-UNAM PAPIME grants PE-105911 and PE-107114, and DGAPA-UNAM PAPIIT grants IN-117214 and IN-117117. LB acknowledges support through a *Cátedra Moshinsky* (2013).
