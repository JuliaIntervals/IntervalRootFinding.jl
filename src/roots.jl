
import IntervalArithmetic: diam, bisect, isnai

export branch_and_prune, Bisection, Newton

struct RootProblem{C, F, G, R, S, T}
    contractor::Type{C}
    f::F
    derivative::G
    region::R
    search_order::Type{S}
    abstol::T
    reltol::T  # TODO
    max_iteration::Int  # TODO
    where_bisect::T
end

RootProblem(f, region ; kwargs...) = RootProblem(f, Root(region, :unkown) ; kwargs...)

function RootProblem(
        f, root::Root ;
        contractor = Newton,
        derivative = nothing,
        search_order = BreadthFirst,
        abstol = 1e-7,
        reltol = NaN,
        max_iteration = 100_000,
        where_bisect = 0.49609375)  # 127//256
    
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
   
function bisect(r::Root, α)
    Y1, Y2 = bisect(root_region(r), α)
    return Root(Y1, :unknown), Root(Y2, :unknown)
end

function process(root_problem, r::Root)
    contracted_root = contract_root(root_problem, r)
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

"""
    roots(f, region ; kwargs...)
"""
function roots(f, region ; kwargs...)
    root_problem = RootProblem(f, region ; kwargs...)
    search = BranchAndPruneSearch(
        root_problem.search_order,
        X -> process(root_problem, X),
        X -> bisect(X, root_problem.where_bisect),
        root_problem.region
    )
    result = bpsearch(search)
    return vcat(result.final_regions, result.unfinished_regions)
end

# TODO Reinstaste support for that
# Acting on complex `Interval`
function _roots(f, Xc::Complex{Interval{T}}, contractor::Type{C},
               search_order::Type{S}, tol::Float64) where {T, C, S <: SearchOrder}

    g = realify(f)
    Y = IntervalBox(reim(Xc)...)
    rts = _roots(g, Root(Y, :unknown), contractor, search_order, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

function _roots(f, Xc::Complex{Interval{T}}, contractor,
               search_order::Type{S}, tol::Float64) where {T, S <: SearchOrder}

    g = realify(f)
    g_prime = x -> ForwardDiff.jacobian(g, x)
    Y = IntervalBox(reim(Xc)...)
    rts = _roots(g, g_prime, Root(Y, :unknown), contractor, search_order, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

function _roots(f, deriv, Xc::Complex{Interval{T}}, contractor,
               search_order::Type{S}, tol::Float64) where {T, S <: SearchOrder}

    g = realify(f)
    g_prime = realify_derivative(deriv)
    Y = IntervalBox(reim(Xc)...)
    rts = _roots(g, g_prime, Root(Y, :unknown), contractor, search_order, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end
