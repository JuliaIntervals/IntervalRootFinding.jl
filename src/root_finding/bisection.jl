
function bisection{T}(f, X::Interval{T}; tolerance=1e-3)

    image = f(X)

    if 0 ∉ image
        return Root{T}[]  # guaranteed that no zero in the interval
    end

    if diam(X) < tolerance
        return [Root{T}(X, :unknown)]
    end

    X1, X2 = bisect(X)

    return [bisection(f, X1, tolerance=tolerance);
            bisection(f, X2, tolerance=tolerance)]

end
