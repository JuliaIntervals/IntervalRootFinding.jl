```@meta
DocTestSetup = quote
    using IntervalRootFinding
end
```
# Decorations and guarantees

When you define a simple function, you will notice that the roots
usually come out with the flag `com` and `NG`.

```jldoctest
julia> roots(x -> x^2 - 0.1, interval(-5, 5))
2-element Vector{Root{Interval{Float64}}}:
 Root([-0.316228, -0.316228]_com_NG, :unique)
 Root([0.316228, 0.316228]_com_NG, :unique)
```

In this case, `com` is the *decoration* of the interval,
and `NG` its *guarantee*.
See the documentation of `IntervalArithmetic.jl` for a detailed description,
here we will go on what they mean in the context of finding roots
of a function.

## Decorations

The decoration tells us whether all operations that affected the interval
were well-behaved.
The two crucial ones for us are `def` (the operations were not continuous)
and `trv` (something horribly wrong happened).

In both cases, it means that the root can not be trusted,
as the hypothesis of the theory we use are not fulfilled.

## Guarantee

The `NG` flags means that non-interval have been mixed with intervals.
Since we can not track the source of the non-interval numbers,
we can not guarantee that they are correct.
For example, `0.1` is famously not parsed as `0.1`,
as `0.1` can not be represented exactly as a binary number
(just like `1/3` can not be represented exactly as a decimal number).

```jldoctest
julia> big(0.1)
0.1000000000000000055511151231257827021181583404541015625
```

If you want the number that you are inputting to be trusted **as is**,
you can use the `@exact` macro from `IntervalArithmetic.jl`,
and the `NG` flag will be avoided in most cases.

```jldoctest exact
julia> @exact f(x) = x^2 - 0.1
f (generic function with 1 method)

julia> roots(f, interval(-5, 5))
2-element Vector{Root{Interval{Float64}}}:
 Root([-0.316228, -0.316228]_com, :unique)
 Root([0.316228, 0.316228]_com, :unique)
```

The `NG` flag can still appear if other computations,
from a library using non-interval "magic" numbers for example,
thus indicating that some non-trusted numbers have been used in the computation.

Moreover, the macro work in such a way that you can still use the defined
function with numbers and get floating point results.

```jldoctest exact
julia> f(0.2)
-0.06
```
## Trust

We are doing our best to really give validated and guaranteed results.
However, in the case that you may make your house explode based on a result
returned by this package,
we would like to remind you that you should not trust the package beyond
the promises of the license:

> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED

All bug reports and pull requests to fix them are however more than welcome.
