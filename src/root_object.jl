"""
    Root

Object representing a possible root inside a given region. The field `status`
is either `:unknown` or `:unique`. If `status` is `:unique` then we know that
there is a unique root of the function in question inside the given region.

Internally the status may also be `:empty` for region guaranteed to contain no
root, however such `Root`s are discarded by default and thus never returned by
the `roots` function.

# Fields
  - `region`: a region (either `Interval` or `SVector` of interval
        representing an interval box) searched for roots.
  - `status`: the status of the region, valid values are `:empty`, `unknown`
        and `:unique`.
"""
struct Root{T}
    region::T
    status::Symbol
end

"""
    root_status(rt)

Return the status of a `Root`.
"""
root_status(rt::Root) = rt.status  # Use root_status since just `status` is too generic

"""
    root_region(rt)

Return the region associated to a `Root`.
"""
root_region(rt::Root) = rt.region

"""
    isunique(rt)

Return whether a `Root` is unique.
"""
isunique(rt::Root{T}) where {T} = (rt.status == :unique)

show(io::IO, rt::Root) = print(io, "Root($(rt.region), :$(rt.status))")

⊆(a::Interval, b::Root) = a ⊆ b.region
⊆(a::Root, b::Root) = a.region ⊆ b.region

IntervalArithmetic.diam(r::Root) = diam_region(root_region(r))
IntervalArithmetic.mag(r::Root) = mag_region(root_region(r))
IntervalArithmetic.mig(r::Root) = mig_region(root_region(r))
IntervalArithmetic.isnai(r::Root) = isnai_region(root_region(r))

function Base.:(==)(r1::Root, r2::Root)
    root_status(r1) != root_status(r2) && return false
    return isequal_region(root_region(r1), root_region(r2))
big(a::Root) = Root(big(a.interval), a.status)

"""
    Base.iterate(r::Root{T})

    Implements the first necessary method for iterate over a Root object.

    Outputs:
    - `r.interval`: interval where some true root of a given function lies.
    - `state`: indicates that the first iterable element of r has been returned.

"""
function Base.iterate(r::Root{T}) where {T}
    state = 1
    return (r.interval, state)
end

"""
    Base.iterate(r::Root{T}, state)

    Implements the second necessary method for iterate over a Root object.

    Inputs:
    - `r`: Root object.
    - `state`: indicates the prior element iterated.

    Outputs:
    - `r.status`: status of the root.
    - `state`: indicates that the second iterable element of r has been returned.
    or 
    - nothing: the iteration has been completed.
"""
function Base.iterate(r::Root{T}, state::Int64) where {T}
    if state == 1
        state = 2
        return (r.status, state)
    end
    # there are no more items to be iterated
    return nothing
end