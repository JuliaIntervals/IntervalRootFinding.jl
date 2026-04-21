
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
    ignored_errors::Vector{DataType}
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

Keyword parameters
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
- `infer_root_type`: When true, use the return type of the function as
    type for the region in the returned roots, avoiding extra conversions
    during the computation.
    Otherwise, use the type of the provided region.
    Default: `true`. Always `false` if the initial region is given as a `Root`.
- `where_bisect`: Value used to bisect the region. It is used to avoid
    bisecting exactly on zero when starting with symmetrical regions,
    often leading to having a solution directly on the boundary of a region,
    which prevent the contractor to prove it's unicity. Default: `127/256`.
-`ignored_errors`: List of exceptions that are ignored during the processing
    of a region. If the error is encoutered, it is discarded and the region is bisected
    further.
    Default: `[IntervalArithmetic.InconclusiveBooleanOperation]`.
"""
function RootProblem(f, region ; infer_root_type = true, kwargs...)
    if infer_root_type
        T = last(InteractiveUtils.@code_typed(f(region)))
        if isconcretetype(T)
            region = convert(T, region)
        else
            @warn "Could not infer the return type of the function (it may be type instable). Got $T"
        end
    end
    RootProblem(f, Root(region, :unkown) ; kwargs...)
end

function RootProblem(
        f, root::Root ;
        contractor = Newton,
        derivative = nothing,
        search_order = BreadthFirst,
        abstol = 1e-7,
        reltol = 0.0,
        max_iteration = 100_000,
        where_bisect = 0.49609375,  # 127//256
        ignored_errors = [IntervalArithmetic.InconclusiveBooleanOperation])
    
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
        convert(Vector{DataType}, ignored_errors)
    )
end

Base.show(io::IO, pb::RootProblem) = print(io, """
    RootProblem
      Contractor: $(pb.contractor)
      Function: $(pb.f)
      Search region: $(root_region(pb.region))
      Search order: $(pb.search_order)
      Absolute tolerance: $(pb.abstol)
      Relative tolerance: $(pb.reltol)
      Maximum iterations: $(pb.max_iteration)
      Ignored errors: $(pb.ignored_errors)"""
)

# TODO Document
# TODO Simplify, only the state is actually needed... Just use it with more methods ?
struct RootSearchState{R <: Root, S <: SearchOrder, B <: BranchAndPruneSearch}
    iteration::Int
    roots::Vector{R}
    state::BranchAndPrune.SearchState{S, R}
    search::B
end

roots(state::RootSearchState) = [leaf.region for leaf in BranchAndPrune.Leaves(state.tree)]
converged_roots(state::RootSearchState) = [leaf.region for leaf in BranchAndPrune.Leaves(state.tree) if leaf.status == :final]
unconverged_roots(state::RootSearchState) = [leaf.region for leaf in state.search_order.working_leaves]

function Base.getproperty(state::RootSearchState, name::Symbol)
    hasfield(RootSearchState, name) && return getfield(state, name)
    hasfield(BranchAndPrune.SearchState, name) && return getfield(state.state, name)
end

function Base.show(io::IO, state::RootSearchState)
    print(io, """RootSearchState
      iteration: $(state.iteration)
      roots:
    """) 
    for rt in state.roots
        println(io, "    $rt")
    end
end

function Base.iterate(root_problem::RootProblem, state = nothing)
    if isnothing(state)
        search = root_search(root_problem)
        bp_iteration = iterate(search)
    else
        search, search_state = state
        bp_iteration = iterate(search, search_state)
    end
    isnothing(bp_iteration) && return nothing
    bp_search_state = bp_iteration[2]
    bp_search_state.iteration > root_problem.max_iteration && return nothing
    state = RootSearchState(
        bp_search_state.iteration,
        [node.region for node in BranchAndPrune.Leaves(bp_search_state.tree)],
        bp_search_state,
        search
    )
    return state, (search, bp_search_state)
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
    contracted = nothing
    try
        contracted = contract(root_problem, root)
    catch err
        !any(isa(err, Err) for Err in root_problem.ignored_errors) && rethrow()
        contracted = Root(root.region, :unknown, :none, (err, backtrace()))
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
            return :store, Root(contracted.region, :unknown, :tolerance, contracted.error)
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
    state = nothing

    for outer state in search end

    result = BranchAndPrune.BranchAndPruneResult(
        state.search_order,
        search.initial_region,
        state.tree
    )

    rts = vcat(result.final_regions, result.unfinished_regions)
    return map(rts) do rt
        if rt.status == :unknown && rt.convergence == :none
            return Root(rt.region, rt.status, :max_iterartion, rt.error)
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
