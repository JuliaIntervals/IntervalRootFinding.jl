
import IntervalArithmetic: diam, isinterior, bisect, isnan

export branch_and_prune, Bisection, Newton
export BreadthFirstSearch, DepthFirstSearch

diam(r::Root) = diam(interval(r))
isnan(X::IntervalBox) = any(isnan.(X))
isnan(r::Root) = isnan(interval(r))

"""
    BreadthFirstSearch <: BreadthFirstBBSearch

Type implementing  the `BreadthFirstBBSearch` interface for interval roots finding.

# Fields:
    - `initial`: region (as a `Root` object) in which roots are searched.
    - `contractor`: contractor to use (`Bisection`, `Newton` or `Krawczyk`)
    - `tol`: tolerance of the search
"""
struct BreadthFirstSearch{R <: Region, C <: Contractor, T <: Real} <: BreadthFirstBBSearch{Root{R}}
    initial::Root{R}
    contractor::C
    tol::T
end

"""
    DepthFirstSearch <: DepthFirstBBSearch

Type implementing the `DepthFirstBBSearch` interface for interval roots finding.

# Fields:
    - `initial`: region (as a `Root` object) in which roots are searched.
    - `contractor`: contractor to use (`Bisection`, `Newton` or `Krawczyk`)
    - `tol`: tolerance of the search
"""
struct DepthFirstSearch{R <: Region, C <: Contractor, T <: Real} <: DepthFirstBBSearch{Root{R}}
    initial::Root{R}
    contractor::C
    tol::T
end

BreadthFirstSearch(X::Region, C::Contractor, tol::Real) = BreadthFirstSearch(Root(X, :unknown), C, tol)
DepthFirstSearch(X::Region, C::Contractor, tol::Real) = DepthFirstSearch(Root(X, :unknown), C, tol)

root_element(search::BBSearch{Root{R}}) where {R <: Region} = search.initial

function bisect(r::Root)
    Y1, Y2 = bisect(interval(r))
    return Root(Y1, :unknown), Root(Y2, :unknown)
end

bisect(::BBSearch, r::Root) = bisect(r::Root)

function process(search::Union{BreadthFirstSearch, DepthFirstSearch}, r::Root)
    contracted_root = search.contractor(r, search.tol)
    status = root_status(contracted_root)

    status == :unique && return :store, contracted_root
    status == :empty && return :discard, contracted_root

    if status == :unknown
        # Avoid infinite division of intervals with singularity
        isnan(contracted_root) && diam(r) < search.tol && return :store, r
        diam(contracted_root) < search.tol && return :store, contracted_root

        return :bisect, r
    else
        error("Unrecognized root status: $status")
    end
end

"""
    branch_and_prune(X, contractor, strategy, tol)

Generic branch and prune routine for finding isolated roots using the given
contractor to determine the status of a given box `X`.

See the documentation of the `roots` function for explanation of the other
arguments.
"""
function branch_and_prune(r::Root, contractor, search, tol)
    iter = search(r, contractor, tol)
    local endstate
    # complete iteration
    for state in iter
        endstate = state
    end
    return data(endstate)
end

const NewtonLike = Union{Type{Newton}, Type{Krawczyk}}
const default_strategy = DepthFirstSearch
const default_tolerance = 1e-15
const default_contractor = Newton

"""
    roots(f, X, contractor=Newton, strategy=BreadthFirstSearch, tol=1e-15)
    roots(f, deriv, X, contractor=Newton, strategy=BreadthFirstSearch, tol=1e-15)
    roots(f, X, contractor, tol)
    roots(f, deriv, X, contractor, tol)

Uses a generic branch and prune routine to find in principle all isolated roots
of a function `f:R^n → R^n` in a region `X`, if the number of roots is finite.

Inputs:
- `f`: function whose roots will be found
- `X`: `Interval` or `IntervalBox` in which roots are searched
- `contractor`: function that, when applied to the function `f`, determines
    the status of a given box `X`. It returns the new box and a symbol indicating
    the status. Current possible values are `Bisection`, `Newton` and `Krawczyk`
- `deriv`: explicit derivative of `f` for `Newton` and `Krawczyk`
- `strategy`: `SearchStrategy` determining the order in which regions are
    processed.
- `tol`: Absolute tolerance. If a region has a diameter smaller than `tol`, it
    is returned with status `:unknown`.

"""
#===
    Default case when `contractor, `strategy` or `tol` is omitted.
===#
function roots(f::Function, X, contractor::Type{C}=default_contractor,
               strategy::Type{S}=default_strategy,
               tol::Float64=default_tolerance) where {C <: Contractor, S <: BBSearch}

    _roots(f, X, contractor, strategy, tol)
end

function roots(f::Function, deriv::Function, X, contractor::Type{C}=default_contractor,
               strategy::Type{S}=default_strategy,
               tol::Float64=default_tolerance) where {C <: Contractor, S <: BBSearch}

    _roots(f, deriv, X, contractor, strategy, tol)
end

function roots(f::Function, X, contractor::Type{C},
               tol::Float64) where {C <: Contractor}

    _roots(f, X, contractor, default_strategy, tol)
end

function roots(f::Function, deriv::Function, X, contractor::Type{C},
               tol::Float64) where {C <: Contractor}

    _roots(f, deriv, X, contractor, default_strategy, tol)
end

#===
    More specific `roots` methods (all parameters are present)
    These functions are called `_roots` to avoid recursive calls.
===#

# For `Bisection` method
function _roots(f, r::Root{T}, ::Type{Bisection},
               strategy::Type{S}, tol::Float64) where {T, S <: BBSearch}

    branch_and_prune(r, Bisection(f), strategy, tol)
end


# For `NewtonLike` acting on `Interval`
function _roots(f, r::Root{Interval{T}}, contractor::NewtonLike,
               strategy::Type{S}, tol::Float64) where {T, S <: BBSearch}

    deriv = x -> ForwardDiff.derivative(f, x)
    _roots(f, deriv, r, contractor, strategy, tol)
end

function _roots(f, deriv, r::Root{Interval{T}}, contractor::NewtonLike,
               strategy::Type{S}, tol::Float64) where {T, S <: BBSearch}

    branch_and_prune(r, contractor(f, deriv), strategy, tol)
end


# For `NewtonLike` acting on `IntervalBox`
function _roots(f, r::Root{IntervalBox{N, T}}, contractor::NewtonLike,
               strategy::Type{S}, tol::Float64) where {N, T, S <: BBSearch}

    deriv = x -> ForwardDiff.jacobian(f, x)
    _roots(f, deriv, r, contractor, strategy, tol)
end

function _roots(f, deriv, r::Root{IntervalBox{N, T}}, contractor::NewtonLike,
               strategy::Type{S}, tol::Float64) where {N, T, S <: BBSearch}

    branch_and_prune(r, contractor(f, deriv), strategy, tol)
end


# Acting on `Interval`
function _roots(f, X::Region, contractor::Type{C},
               strategy::Type{S}, tol::Float64) where {C <: Contractor, S <: BBSearch}

    _roots(f, Root(X, :unknown), contractor, strategy, tol)
end

function _roots(f, deriv, X::Region, contractor::Type{C},
               strategy::Type{S}, tol::Float64) where {C <: Contractor, S <: BBSearch}

    _roots(f, deriv, Root(X, :unknown), contractor, strategy, tol)
end


# Acting on `Vector` of `Root`
function _roots(f, V::Vector{Root{T}}, contractor::Type{C},
               strategy::Type{S}, tol::Float64) where {T, C <: Contractor, S <: BBSearch}

    vcat(_roots.(f, V, contractor, strategy, tol)...)
end

function _roots(f, deriv, V::Vector{Root{T}}, contractor::Type{C},
               strategy::Type{S}, tol::Float64) where {T, C <: Contractor, S <: BBSearch}

    vcat(_roots.(f, deriv, V, contractor, strategy, tol)...)
end


# Acting on complex `Interval`
function _roots(f, Xc::Complex{Interval{T}}, contractor::Type{C},
               strategy::Type{S}, tol::Float64) where {T, C <: Contractor, S <: BBSearch}

    g = realify(f)
    Y = IntervalBox(reim(Xc)...)
    rts = _roots(g, Root(Y, :unknown), contractor, strategy, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

function _roots(f, Xc::Complex{Interval{T}}, contractor::NewtonLike,
               strategy::Type{S}, tol::Float64) where {T, S <: BBSearch}

    g = realify(f)
    g_prime = x -> ForwardDiff.jacobian(g, x)
    Y = IntervalBox(reim(Xc)...)
    rts = _roots(g, g_prime, Root(Y, :unknown), contractor, strategy, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

function _roots(f, deriv, Xc::Complex{Interval{T}}, contractor::NewtonLike,
               strategy::Type{S}, tol::Float64) where {T, S <: BBSearch}

    g = realify(f)
    g_prime = realify_derivative(deriv)
    Y = IntervalBox(reim(Xc)...)
    rts = _roots(g, g_prime, Root(Y, :unknown), contractor, strategy, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end
