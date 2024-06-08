# IntervalRootFinding.jl

[![Build Status](https://github.com/JuliaIntervals/IntervalRootFinding.jl/workflows/CI/badge.svg)](https://github.com/JuliaIntervals/IntervalRootFinding.jl/actions/workflows/CI.yml)

[![coverage](https://codecov.io/gh/JuliaIntervals/IntervalRootFinding.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaIntervals/IntervalRootFinding.jl)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaintervals.github.io/IntervalRootFinding.jl/latest/)

This package provides guaranteed methods for finding **roots** of functions, i.e. solutions to the equation `f(x) == 0` for a function `f`.
To do so, it uses methods from interval analysis, using interval arithmetic from the [`IntervalArithmetic.jl`](https://github.com/JuliaIntervals/IntervalArithmetic.jl) package by the same authors.

## Basic usage examples

The basic function is `roots`. A standard Julia function and an interval is provided and the `roots` function return a list of intervals containing *all* roots of the function located in the starting interval.

```jl
julia> using IntervalArithmetic, IntervalArithmetic.Symbols, IntervalRootFinding

julia> f(x) = sin(x) - 0.1*x^2 + 1
f (generic function with 1 method)

julia> roots(f, -10..10)
4-element Vector{Root{Interval{Float64}}}:
 Root([-4.42654, -4.42653]_com_NG, :unique)
 Root([-3.10682, -3.10681]_com_NG, :unique)
 Root([-1.08205, -1.08204]_com_NG, :unique)
 Root([3.14959, 3.1496]_com_NG, :unique)
```

The `:unique` status tell us, in addition, that each listed region contains exactly one root. The other possible status is `:unknown`, which corresponds to intervals that may contain zero, one, or more roots - no guarantee is provided for these intervals.

These results are represented in the following plot, the region containing roots being in green. The inset show a close-up of one of the roots:

![basic usage](docs/src/basic_usage.png)

The full documentation is available [here](https://juliaintervals.github.io/IntervalRootFinding.jl/latest/).


### Note about guarantee

In the example, the interval of the roots have the `NG` flag, meaning that they are not strictly guaranteed to be correct.
This happens when numbers and interval are mixed,
because IntervalArithmetic does not know where the numbers come from,
and can therefore not guarantee that they are properly rounded.

The NG flag can usually be avoided by taking the responsibility
to provide exactly the right number using the `@exact` macro.

```julia
julia> @exact g(x) = sin(x) - 0.1*x^2 + 1
g (generic function with 1 method)

julia> roots(g, -10..10)
4-element Vector{Root{Interval{Float64}}}:
 Root([-4.42654, -4.42653]_com, :unique)
 Root([-3.10682, -3.10681]_com, :unique)
 Root([-1.08205, -1.08204]_com, :unique)
 Root([3.14959, 3.1496]_com, :unique)
 ```

More details about the interval decoration system and what it means for root finding can be found in the documentation.

## Authors
- [Luis Benet](http://www.cicc.unam.mx/~benet/), Instituto de Ciencias Físicas,
Universidad Nacional Autónoma de México (UNAM)
- [David P. Sanders](http://sistemas.fciencias.unam.mx/~dsanders),
Departamento de Física, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM)

## Acknowledgements ##

Financial support is acknowledged from DGAPA-UNAM PAPIME grants PE-105911 and PE-107114, and DGAPA-UNAM PAPIIT grants IN-117214 and IN-117117. LB acknowledges support through a *Cátedra Moshinsky* (2013).
