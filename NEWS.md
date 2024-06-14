# Updates to `IntervalRootFinding.jl`

## v0.6

Compatibility with IntervalArithmetic v0.22, which fundamentally changes the package.
- Decorated intervals must be used (which is the default `Interval` in the latest releases of IntervalArithmetic)
- The signature of the `roots` function changed to `roots(f::Function, X::Union{Interval, AbstractVector} ; kwargs...)`, with the following consequences
    - Manual derivatives and contractors must now always be passed as keyword arguments. This greatly simplify the internal logic of the `roots` function.
    - No more `IntervalBox`. Multidimensional problem are specified by returning a vector of intervals, and giving a vector of intervals as initial search region.
    - Normal vectors of intervals are supported. They are significantly (3x) time slower for a simple 2D problem than SVector, but more convenient.
- Features of the packages outside the `roots` function are unsupported for the time being (e.g. `Slopes` and quadratic equations). If you need them, please open an issue.
- `roots` works with `@exact`.


## v0.2

The package has been completely rewritten for v0.2.
It now contains basic bisection and Newton interval routines for functions with an arbitrary number of parameters.
