
import IntervalArithmetic: diam, bisect, isnan

export branch_and_prune, Bisection, Newton

diam(r::Root) = diam(interval(r))
isnan(r::Root) = isnan(interval(r))

struct RootProblem{T}
    abstol::T
end
   
function bisect(r::Root)
    Y1, Y2 = bisect(interval(r))
    return Root(Y1, :unknown), Root(Y2, :unknown)
end

function process(contractor, root_problem, r::Root)
    contracted_root = contract(contractor, r)
    refined_root = refine(contractor, contracted_root, root_problem)

    status = root_status(refined_root)

    status == :unique && return :store, refined_root
    status == :empty && return :prune, refined_root

    if status == :unknown
        # Avoid infinite division of intervals with singularity
        isnan(refined_root) && diam(r) < root_problem.abstol && return :store, r
        diam(refined_root) < root_problem.abstol && return :store, refined_root

        return :branch, r
    else
        error("Unrecognized root status: $status")
    end
end

"""
    branch_and_prune(X, contractor, search_order, tol)

Generic branch and prune routine for finding isolated roots using the given
contractor to determine the status of a given box `X`.

See the documentation of the `roots` function for explanation of the other
arguments.
"""
function branch_and_prune(r::Root, contractor, search_order, tol)
    root_problem = RootProblem(tol)
    search = BranchAndPruneSearch(
        search_order,
        X -> process(contractor, root_problem, X),
        bisect,
        r
    )
    result = bpsearch(search)
    return vcat(result.final_regions, result.unfinished_regions)
end

const NewtonLike = Union{Type{Newton}, Type{Krawczyk}}
const default_search_order = DepthFirst
const default_tolerance = 1e-7
const default_contractor = Newton

"""
    roots(f, X, contractor=Newton, search_order=BreadthFirst, tol=1e-15)
    roots(f, deriv, X, contractor=Newton, search_order=BreadthFirst, tol=1e-15)
    roots(f, X, contractor, tol)
    roots(f, deriv, X, contractor, tol)

Uses a generic branch and prune routine to find in principle all isolated roots
of a function `f:R^n → R^n` in a region `X`, if the number of roots is finite.

Inputs:
  - `f`: function whose roots will be found
  - `X`: `Interval` or `SVector` of `Interval` in which roots are searched
  - `contractor`: function that, when applied to the function `f`, determines
    the status of a given box `X`. It returns the new box and a symbol
    indicating the status. Current possible values are `Bisection`, `Newton`
    and `Krawczyk`
  - `deriv`: explicit derivative of `f` for `Newton` and `Krawczyk`
  - `search_order`: `SearchStrategy` determining the order in which regions are
    processed.
  - `tol`: Absolute tolerance. If a region has a diameter smaller than `tol`,
    it is returned with status `:unknown`.

"""
function roots(f::Function, X, contractor::Type{C}=default_contractor,
               search_order::Type{S}=default_search_order,
               tol::Float64=default_tolerance) where {C <: AbstractContractor, S <: SearchOrder}

    _roots(f, X, contractor, search_order, tol)
end

function roots(f::Function, deriv::Function, X, contractor::Type{C}=default_contractor,
               search_order::Type{S}=default_search_order,
               tol::Float64=default_tolerance) where {C <: AbstractContractor, S <: SearchOrder}

    _roots(f, deriv, X, contractor, search_order, tol)
end

function roots(f::Function, X, contractor::Type{C},
               tol::Float64) where {C <: AbstractContractor}

    _roots(f, X, contractor, default_search_order, tol)
end

function roots(f::Function, deriv::Function, X, contractor::Type{C},
               tol::Float64) where {C <: AbstractContractor}

    _roots(f, deriv, X, contractor, default_search_order, tol)
end

#===
    More specific `roots` methods (all parameters are present)
    These functions are called `_roots` to avoid recursive calls.
===#

# For `Bisection` method
function _roots(f, r::Root{T}, ::Type{Bisection},
               search_order::Type{S}, tol::Float64) where {T, S <: SearchOrder}

    branch_and_prune(r, Bisection(f), search_order, tol)
end


# For `NewtonLike` acting on `Interval`
function _roots(f, r::Root{Interval{T}}, contractor::NewtonLike,
               search_order::Type{S}, tol::Float64) where {T, S <: SearchOrder}

    deriv = x -> ForwardDiff.derivative(f, x)
    _roots(f, deriv, r, contractor, search_order, tol)
end

function _roots(f, deriv, r::Root{Interval{T}}, contractor::NewtonLike,
               search_order::Type{S}, tol::Float64) where {T, S <: SearchOrder}

    branch_and_prune(r, contractor(f, deriv), search_order, tol)
end


# For `NewtonLike` acting on `IntervalBox`
function _roots(f, r::Root{<:SVector}, contractor::NewtonLike,
               search_order::Type{<:SearchOrder}, tol::Float64)

    deriv = x -> ForwardDiff.jacobian(f, x)
    _roots(f, deriv, r, contractor, search_order, tol)
end

function _roots(f, deriv, r::Root{<:SVector}, contractor::NewtonLike,
               search_order::Type{<:SearchOrder}, tol::Float64)

    branch_and_prune(r, contractor(f, deriv), search_order, tol)
end


# Acting on `Interval`
function _roots(f, X::Region, contractor::Type{C},
               search_order::Type{S}, tol::Float64) where {C <: AbstractContractor, S <: SearchOrder}

    _roots(f, Root(X, :unknown), contractor, search_order, tol)
end

function _roots(f, deriv, X::Region, contractor::Type{C},
               search_order::Type{S}, tol::Float64) where {C <: AbstractContractor, S <: SearchOrder}

    _roots(f, deriv, Root(X, :unknown), contractor, search_order, tol)
end


# Acting on `Vector` of `Root`
function _roots(f, V::Vector{Root{T}}, contractor::Type{C},
               search_order::Type{S}, tol::Float64) where {T, C <: AbstractContractor, S <: SearchOrder}

    vcat(_roots.(f, V, contractor, search_order, tol)...)
end

function _roots(f, deriv, V::Vector{Root{T}}, contractor::Type{C},
               search_order::Type{S}, tol::Float64) where {T, C <: AbstractContractor, S <: SearchOrder}

    vcat(_roots.(f, deriv, V, contractor, search_order, tol)...)
end


# Acting on complex `Interval`
function _roots(f, Xc::Complex{Interval{T}}, contractor::Type{C},
               search_order::Type{S}, tol::Float64) where {T, C <: AbstractContractor, S <: SearchOrder}

    g = realify(f)
    Y = IntervalBox(reim(Xc)...)
    rts = _roots(g, Root(Y, :unknown), contractor, search_order, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

function _roots(f, Xc::Complex{Interval{T}}, contractor::NewtonLike,
               search_order::Type{S}, tol::Float64) where {T, S <: SearchOrder}

    g = realify(f)
    g_prime = x -> ForwardDiff.jacobian(g, x)
    Y = IntervalBox(reim(Xc)...)
    rts = _roots(g, g_prime, Root(Y, :unknown), contractor, search_order, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

function _roots(f, deriv, Xc::Complex{Interval{T}}, contractor::NewtonLike,
               search_order::Type{S}, tol::Float64) where {T, S <: SearchOrder}

    g = realify(f)
    g_prime = realify_derivative(deriv)
    Y = IntervalBox(reim(Xc)...)
    rts = _roots(g, g_prime, Root(Y, :unknown), contractor, search_order, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end
