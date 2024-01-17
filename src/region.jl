function intersect_region(X::Interval, Y::Interval)
    intersection = intersect_interval(bareinterval(X), bareinterval(Y))
    dec = min(decoration(X), decoration(Y))
    guarantee = isguaranteed(X) && isguaranteed(Y)
    return IntervalArithmetic._unsafe_interval(intersection, dec, guarantee)
end

intersect_region(X::AbstractVector, Y::AbstractVector) = intersect_region.(X, Y)

isempty_region(X::Interval) = isempty_interval(X)
isempty_region(X::AbstractVector) = any(isempty_region.(X))

isequal_region(X::Interval, Y::Interval) = isequal_interval(X, Y)
isequal_region(X::AbstractVector, Y::AbstractVector) = all(isequal_region.(X, Y))

isbounded_region(X::Interval) = isbounded(X)
isbounded_region(X::AbstractVector) = all(isbounded.(X))

isnai_region(X::Interval) = isnai(X)
isnai_region(X::AbstractVector) = any(isnai.(X))

diam_region(X::Interval) = diam(X)
diam_region(X::AbstractVector) = maximum(diam.(X))

function bisect(X::Interval, α)
    m = mid(X, α)
    return (interval(inf(X), m), interval(m, sup(X)))
end

function bisect(X::AbstractVector, α)
    X1 = copy(X)
    X2 = copy(X)

    i = argmax(diam.(X))
    x1, x2 = bisect(X[i], α)
    X1[i] = x1
    X2[i] = x2
    return X1, X2
end