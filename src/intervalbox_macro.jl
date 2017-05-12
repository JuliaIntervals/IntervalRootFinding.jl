using ValidatedNumerics


doc"""
Macro to make a version of a multidimensional function that acts on
an `IntervalBox` and returns an `IntervalBox`.

The original n-argument function, the `IntervalBox` version, and a version
taking and returning an array, for use with automatic differentiation via `ForwardDiff`, are created.

Example:

```
julia v0.5> @intervalbox f(x, y) = (x + y, x - y)
f (generic function with 3 methods)

julia v0.5> X = IntervalBox(1..1, 2..2)
[1, 1] × [2, 2]

julia v0.5> f(X)
[3, 3] × [-1, -1]

julia v0.5> ForwardDiff.jacobian(f, [1, 2])
2×2 Array{Int64,2}:
 1   1
 1  -1
```

(No significant error checking is performed.)
"""

macro intervalbox(ex)
    # @show ex
    # @show ex.head
    if ex.head != :(=)
        throw(ArgumentError("Not a function definition."))
    end

    if ex.args[1].head != :call
        throw(ArgumentError("Not a function definition."))
    end

    f = ex.args[1].args[1]  # function name
    f = esc(f)
    ex.args[1].args[1] = f  # escape the function name in the original code


    new_ex = Expr(:block)

    push!(new_ex.args, ex)  # the function call
    push!(new_ex.args,  :( ($f)(X::Vector) = [$f(X...)...] ) )
    push!(new_ex.args,  :( ($f)(X::IntervalBox) = IntervalBox( $f(X...)...) ) )

    #@show new_ex
    new_ex
end
