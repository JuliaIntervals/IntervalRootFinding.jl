struct Region{T, N}
    intervals::T
end

Region(X::T) where {T <: Interval} = Region{T, 1}(X)

function Region(Xs::AbstractVector{T}) where {T <: Interval}
    N = length(Xs)
    N == 1 && return Region{T, 1}(only(Xs))
    return Region{T, N}(Xs) 
end

function Base.intersect(X::Region{<:Any, 1}, Y::Region{<:Any, 1})
    x = only(X.intervals)
    y = only(Y.intervals)
    intersection = intersect_interval(bareinterval(x), bareinterval(x))
    dec = min(decoration(x), decoration(y))
    guarantee = isguaranteed(x) && isguaranteed(y)
    decorated = IntervalArithmetic._unsafe_interval(intersection, dec, guarantee)
    return Region(decorated)
end

function Base.intersect(X::Region{<:Any, N}, Y::Region{<:Any, N}) where N
    return Region(intersect.(X.intervals, Y.intervals))
end
