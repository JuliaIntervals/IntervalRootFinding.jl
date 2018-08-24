"""
Object representing a possible root inside a given `interval`.
`status` is either `:unknown` or `:unique`.
If `status` is `:unique` then we know that there is a unique root of the
function in question inside the given `interval`.
"""
struct Root{T}
    interval::T
    status::Symbol
end

interval(x::Root) = x.interval

# use root_status since just `status` is too generic
root_status(x::Root) = x.status

show(io::IO, root::Root) = print(io, "Root($(root.interval), :$(root.status))")

isunique(root::Root{T}) where {T} = (root.status == :unique)

⊆(a::Interval, b::Root) = a ⊆ b.interval   # the Root object has the interval in the first entry
⊆(a::Root, b::Root) = a.interval ⊆ b.interval

big(a::Root) = Root(big(a.interval), a.status)
