"""
    Root

Object representing a possible root inside a given region. The field `status`
is either `:unknown` or `:unique`. If `status` is `:unique` then we know that
there is a unique root of the function in question inside the given region.

Internally the status may also be `:empty` for region guaranteed to contain no
root, however such `Root`s are discarded by default and thus never returned by
the `roots` function.

# Fields
  - `interval`: a region (either `Interval` or `SVector` of interval
        representing an interval box) searched for roots.
  - `status`: the status of the region, valid values are `:empty`, `unknown`
        and `:unique`.
"""
struct Root{T, N}
    region::IntervalRegion{T, N}
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

diam(r::Root) = diam(interval(r))
isnai(r::Root) = isnai(interval(r))