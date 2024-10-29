
import IntervalArithmetic: diam, bisect, isnai

export branch_and_prune, Bisection, Newton

struct RootProblem{C, F, G, R, S, T}
    contractor::Type{C}
    f::F
    derivative::G
    region::R
    search_order::Type{S}
    abstol::T
    reltol::T
    max_iteration::Int
    where_bisect::T
end

"""
    RootProblem(f::Function, search_region ; kwargs...)

Setup a `RootProblem` for searching the roots (also known as zeros)
of the function `f` in the given search region.

`search_region` must be either an `Interval` if `f` is scalar
or a vector of `Interval`s if `f` is vector-valued.

The returned `RootProblem` is an iterator that give access to the internal
state of the search during the iteration,
allowing to add callbacks and logging to the search.

Parameters
==========
- `contractor`: Contractor used to determine the status of a region.
    Must be either `Newton`, `Krawczyk`, or `Bisection`. `Bisection` do not require
    to compute the derivative of the function, but can also never guarantee the
    existence of a root. Default: `Newton`.
- `derivative`: Explicit derivative of the function (or its jacobian for
    vector-valued functions) used by the `Newton` and `Krawczyk` contractors.
    Default: `nothing` (the derivative is computed automatically using ForwardDiff.jl).
- `search_order`: Order in which the sub-regions are searched.
    `BreadthFirst` (visit the largest regions first) and `DepthFirst`
    (visit the smallest regions first) are supported. Default: `BreadthFirst`.
- `abstol`: Absolute tolerance. The search is stopped when all dimension
    of the remaining regions are smaller than `abstol`. Default: 1e-7.
- `where_bisect`: Value used to bisect the region. It is used to avoid
    bisecting exactly on zero when starting with symmetrical regions,
    often leading to having a solution directly on the boundary of a region,
    which prevent the contractor to prove it's unicity. Default: 127/256.

`reltol` and `max_iteration` are currently ignored.
"""
RootProblem(f, region ; kwargs...) = RootProblem(f, Root(region, :unkown) ; kwargs...)

function RootProblem(
        f, root::Root ;
        contractor = Newton,
        derivative = nothing,
        search_order = BreadthFirst,
        abstol = 1e-7,
        reltol = nothing,
        max_iteration = nothing,
        where_bisect = 0.49609375)  # 127//256
    
    if !isnothing(reltol) || !isnothing(max_iteration)
        throw(
            ArgumentError("reltol and max_iteration not yet implemented")
        )
    else
        reltol = 0.0
        max_iteration = -1
    end
    
    N = length(root_region(root))
    if isnothing(derivative)
        if N == 1
            derivative = x -> ForwardDiff.derivative(f, x)
        else
            derivative = x -> ForwardDiff.jacobian(f, x)
        end
    end

    RootProblem(
        contractor,
        f,
        derivative,
        root,
        search_order,
        abstol,
        reltol,
        max_iteration,
        where_bisect
    )
end

function Base.iterate(root_problem::RootProblem, state = nothing)
    if isnothing(state)
        search = root_search(root_problem)
        iteration = iterate(search)
    else
        search, search_state = state
        iteration = iterate(search, search_state)
    end
    isnothing(iteration) && return nothing
    value, new_search_state = iteration
    return value, (search, new_search_state)
end
   
function bisect_region(r::Root, α)
    Y1, Y2 = bisect_region(root_region(r), α)
    return Root(Y1, :unknown), Root(Y2, :unknown)
end

function process(root_problem, r::Root)
    contracted_root = contract(root_problem, r)
    refined_root = refine(root_problem, contracted_root)

    status = root_status(refined_root)

    status == :unique && return :store, refined_root
    status == :empty && return :prune, refined_root

    if status == :unknown
        # Avoid infinite division of intervals with singularity
        isnai(refined_root) && diam(r) < root_problem.abstol && return :store, r
        diam(refined_root) < root_problem.abstol && return :store, refined_root

        return :branch, r
    else
        error("Unrecognized root status: $status")
    end
end

root_search(root_problem::RootProblem) = BranchAndPruneSearch(
    root_problem.search_order,
    X -> process(root_problem, X),
    X -> bisect_region(X, root_problem.where_bisect),
    root_problem.region
)

"""
    roots(f::Function, search_region ; kwargs...)

Return the roots (also known as zeros) of the function `f`
contained in the given search region
(either an `Interval` for a scalar function or vector of `Interval`s for a
vector valued function),
together with a status.

The status of the returned regions can be either of
- `:unique`: the region provably contains exactly one root of `f`
- `:unknown`: the region may contain any number of roots (potentially none)

The parts of the search region that are not covered by any of the returned
roots are guaranteed to contain no root of the function.

For information about the optional search parameters,
see [`RootProblem`](@ref).
"""
function roots(f, region ; kwargs...)
    search = root_search(RootProblem(f, region ; kwargs...))
    result = bpsearch(search)
    return vcat(result.final_regions, result.unfinished_regions)
end

# Acting on complex `Interval`
function roots(f, region::Complex{<:Interval} ; derivative = nothing, contractor = Newton, kwargs...)
    g = realify(f)
    X = [real(region), imag(region)]

    contractor != Bisection && isnothing(derivative) && throw(
        ArgumentError("when using complex numbers, the derivative must be " *
        "given explicitly as ForwardDiff does not support complex numbers.")
    )

    dg = realify_derivative(derivative)

    rts = roots(g, X ; derivative = dg, contractor, kwargs...)
    return [Root(Complex(root_region(rt)...), root_status(rt)) for rt in rts]
end
