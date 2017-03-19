
doc"""
    bisection(f, X; tolerance=1e-3)

Find possible roots of the function `f` inside the `Interval` or `IntervalBox` `X`.
"""
function bisection(f, X::Union{Interval, IntervalBox}; tolerance=1e-3)

    image = f(X)

    if !(zero(X) ⊆ image)
        return Root{T}[]
    end

    if diam(X) < tolerance
        return [Root{T}(X, :unknown)]
    end

    X1, X2 = bisect(X)

    return [bisection(f, X1, tolerance=tolerance);
            bisection(f, X2, tolerance=tolerance)]

end
