
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
    bisect_on_error::Bool
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
- `abstol`: Absolute tolerance. The search is stopped when all dimensions
    of the remaining regions are smaller than `abstol`. Default: `1e-7`.
- `reltol`: Relative tolerance. The search is stopped when all dimensions
    of the remaining regions are smaller than `reltol * mag(interval)`.
    Default: `0.0`.
- `max_iteration`: The maximum number of iteration, which also corresponds to
    the maximum number of bisections allowed. Default: `100_000`.
- `where_bisect`: Value used to bisect the region. It is used to avoid
    bisecting exactly on zero when starting with symmetrical regions,
    often leading to having a solution directly on the boundary of a region,
    which prevent the contractor to prove it's unicity. Default: `127/256`.
-`bisect_on_error`: Whether a region that errors when the function is applied to
    to it should be bisected. If false, when an error happen, the root search is
    interrupted. Default: true.
"""
RootProblem(f, region ; kwargs...) = RootProblem(f, Root(region, :unkown) ; kwargs...)

function RootProblem(
        f, root::Root ;
        contractor = Newton,
        derivative = nothing,
        search_order = BreadthFirst,
        abstol = 1e-7,
        reltol = 0.0,
        max_iteration = 100_000,
        where_bisect = 0.49609375,  # 127//256
        bisect_on_error = true)
    
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
        where_bisect,
        bisect_on_error
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

function under_tolerance(root_problem, root::Root)
    d = diam(root)
    d < root_problem.abstol && return true
    d / mag(root) < root_problem.reltol && return true
    return false
end

function process(root_problem, root::Root)
    if root_problem.bisect_on_error
        try
            contracted = contract(root_problem, root)
        catch
            contracted = Root(root.region, :unknown, :none, true)
        end
    else
        contracted = contract(root_problem, root)
    end

    status = root_status(contracted)

    if status == :unique
        refined_root = refine(root_problem, contracted)
        return :store, refined_root
    end

    status == :empty && return :prune, root

    if status == :unknown
        # Avoid infinite division of intervals with singularity
        if istrivial(contracted.region) && under_tolerance(root_problem, root)
            return :store, Root(root.region, :unknown, :tolerance)
        end

        if under_tolerance(root_problem, contracted)
            return :store, Root(contracted.region, :unknown, :tolerance, contracted.errored)
        end
        
        return :branch, root
    else
        throw(ArgumentError("unrecognized root status: $status"))
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
    problem = RootProblem(f, region ; kwargs...)
    search = root_search(problem)
    endstate = nothing

    for (iter, state) in enumerate(search)
        endstate = state
        iter >= problem.max_iteration && break
    end

    result = BranchAndPrune.BranchAndPruneResult(
        endstate.search_order,
        search.initial_region,
        endstate.tree
    )

    rts = vcat(result.final_regions, result.unfinished_regions)
    return map(rts) do rt
        if rt.status == :unknown && rt.convergence == :none
            return Root(rt.region, rt.status, :max_iter, rt.errored)
        end
        return rt
    end
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
