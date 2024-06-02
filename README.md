<h1 align="center">
IntervalRootFinding

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaintervals.github.io/IntervalRootFinding.jl/stable)
[![Build Status](https://github.com/JuliaIntervals/IntervalRootFinding.jl/workflows/CI/badge.svg)](https://github.com/JuliaIntervals/IntervalRootFinding.jl/actions/workflows/CI.yml)
[![coverage](https://codecov.io/gh/JuliaIntervals/IntervalRootFinding.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaIntervals/IntervalRootFinding.jl)
</h1>

IntervalRootFinding.jl is a Julia package for finding the roots of functions, i.e. solutions to the equation $f(x) = 0$.
To do so, it uses interval arithmetic from the [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl) library.

## Documentation

The official documentation is available online: https://juliaintervals.github.io/IntervalRootFinding.jl/stable.

## Installation

The IntervalRootFinding.jl package requires to [install Julia](https://julialang.org/downloads/).

Then, start Julia and execute the following command in the REPL:

```julia
using Pkg; Pkg.add("IntervalRootFinding")
```

## Basic usage examples

The basic function is `roots`. Given a standard Julia function and an interval, the `roots` function returns a list of intervals containing all roots of the function located in the prescribed interval.

```julia
julia> using IntervalArithmetic, IntervalRootFinding

julia> f(x) = sin(x) - 0.1*x^2 + 1
f (generic function with 1 method)

julia> roots(f, -10..10)
4-element Array{Root{Interval{Float64}},1}:
 Root([3.14959, 3.1496], :unique)
 Root([-4.42654, -4.42653], :unique)
 Root([-1.08205, -1.08204], :unique)
 Root([-3.10682, -3.10681], :unique)
```

The `:unique` status indicates that each listed interval contains exactly one root. The other possible status is `:unknown`, which corresponds to intervals that may contain zero, one, or more roots (no guarantee is provided for these intervals).

These results are represented in the following plot, the region containing roots being in green.

![basic usage](docs/src/basic_usage.png)
