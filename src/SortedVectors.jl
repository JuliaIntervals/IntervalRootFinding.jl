module SortedVectors

import Base: push!, insert!, deleteat!, resize!, getindex, show, length

export SortedVector

struct SortedVector{T,F}
    data::Vector{T}
    by::F

    function SortedVector(data::Vector{T}, by::F) where {T,F}
        new{T,F}(sort(data), by)
    end
end

SortedVector(data::Vector) = SortedVector(data, identity)

function show(io::IO, v::SortedVector)
    print(io, "SortedVector($(v.data))")
end



function push!{T}(v::SortedVector{T}, x::T)
    push!(v.data, x)
    sort!(v.data, by=v.by)
    return v
end

getindex(v::SortedVector, i::Int) = v.data[i]
length(v::SortedVector) = length(v.data)


testdata = SortedVector([(3, 1), (2, 2)], x->x[1])

testdata

push!(testdata, (1, 3))


function insert!{T}(v::SortedVector{T}, i::Int, x::T)
    push!(v.data, i, x)
    return v
end

function insert!{T}(v::SortedVector{T}, x::T)
    i = binary_search(v, x)
    insert!(v.data, i, x)
    return v
end

function deleteat!(v::SortedVector, i::Int)
    deleteat!(v.data, i)
    return v
end

function resize!(v::SortedVector, n::Int)
    resize!(v.data, n)
    return v
end


"""
Do binary search for item `x` in *sorted* vector `v`.
Returns the lower bound for the position of `x` in `v`.
"""
function binary_search(v, x)
    a, b = 1, length(v)

    if x < v[a]
        return 1

    elseif x > v[b]
        return b
    end

    m = (a + b) ÷ 2  # mid-point

    while abs(a - b) > 1

        @show a, b, v[a], v[b]

        if v[m] == x
            return m
        end

        if x < v[m]
            b = m
        else
            a = m
        end
    end

    @show a, b, v[a], v[b]

    if v[a] == x
        return a
    end

    if v[b] == x
        return b
    end

    return b  # a is a lower_bound; insert in position b
end

v = SortedVector([1, 3, 6, 7, 9, 10])
# binary_search(v, 3)
# binary_search(v, 4)
# binary_search(v, 1)

v2 = SortedVector([2.1, 2.5, 3.0, 4.0])

end
