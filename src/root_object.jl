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
  - `convergence`: the convergence status of the region. It is always `:converged`
        for roots with status `:unique`,
        and can be either `:max_iter` or `:tolerance` for roots with status `:unknown`,
        depending on whether they stopped being processing due to reaching
        the maximum number of iteration or the tolerance, respectively.
  - `errored`: whether an error was encounter during the processing of this region.
        Errors can be raised explicitly when encountered by setting
        `bisect_on_error` to `flase` in the RootProblem.
"""
struct Root{T}
    region::T
    status::Symbol
    convergence::Symbol
    errored::Bool
end

function Root(region, status::Symbol, convergence = :none, errored = false)
    return Root(region, status, convergence, errored)
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

function show(io::IO, rt::Root)
    print(io, "Root($(rt.region), :$(rt.status))")
    if rt.status == :unknown
        if rt.convergence == :tolerance
            print(io, "\n    Not converged: region size smaller than the tolerance")
        elseif rt.convergence == :max_iter 
            print(io, "\n    Not converged: reached maximal number of iterations")
        else
            print(io, "\n    Not converged: unknown reason $(rt.convergence)")
        end

        if rt.errored
            print(io, "\n    Warning: an error was encountered in during computation")
        end
    end
end

⊆(a::Interval, b::Root) = a ⊆ b.region
⊆(a::Root, b::Root) = a.region ⊆ b.region

IntervalArithmetic.diam(r::Root) = diam_region(root_region(r))
IntervalArithmetic.mag(r::Root) = mag_region(root_region(r))
IntervalArithmetic.mig(r::Root) = mig_region(root_region(r))
IntervalArithmetic.isnai(r::Root) = isnai_region(root_region(r))

function Base.:(==)(r1::Root, r2::Root)
    root_status(r1) != root_status(r2) && return false
    return isequal_region(root_region(r1), root_region(r2))
end

big(a::Root) = Root(big(a.interval), a.status)

"""
    Base.iterate(r::Root [, state])

Return successively the root region and status,
allowing to unpack the root object as `region, status = root`.
"""
function Base.iterate(r::Root{T}, state::Integer=1) where {T}
    state == 1 && return (r.region, 2)
    state == 2 && return (r.status, 3)
    return nothing
end
