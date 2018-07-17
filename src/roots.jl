
import IntervalArithmetic: diam, isinterior

export branch_and_prune, Bisection, Newton

diam(x::Root) = diam(x.interval)

"""
    branch_and_prune(X, contractor, strategy, tol)

Generic branch and prune routine for finding isolated roots using the `contract`
function as the contractor. The argument `contractor` is function that
determines the status of a given box `X`. It returns a new contracted box and
a symbol indicating its status.

See the documentation of the `roots` function for explanation of the other
arguments.
"""
function branch_and_prune(X, contractor, strategy, tol)
    iter = RootSearch(X, contractor, strategy, tol)
    local state
    # complete iteration
    for state in iter end
    return state.outputs
end

export recursively_branch_and_prune

function recursively_branch_and_prune(h, X, contractor=BisectionContractor, final_tol=1e-14)
    tol = 2
    roots = branch_and_prune(h, X, IntervalRootFinding.BisectionContractor, tol)

    while tol > 1e-14
       tol /= 2
       roots = branch_and_prune(h, roots, IntervalRootFinding.BisectionContractor, tol)
    end

    return roots
end


const IntervalLike{T} = Union{Interval{T}, IntervalBox{T}}
const NewtonLike = Union{Type{Newton}, Type{Krawczyk}}

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
    is returned with status `:unkown`.

"""

const default_strategy = DepthFirstSearch
const default_tolerance = 1e-15
const default_contractor = Newton

#===
    Default case when `contractor, `strategy` or `tol` is omitted.
===#
function roots(f::Function, X, contractor::Type{C}=default_contractor,
               strategy::SearchStrategy=default_strategy,
               tol::Float64=default_tolerance) where {C <: Contractor}

    _roots(f, X, contractor, strategy, tol)
end

function roots(f::Function, deriv::Function, X, contractor::Type{C}=default_contractor,
               strategy::SearchStrategy=default_strategy,
               tol::Float64=default_tolerance) where {C <: Contractor}

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
function _roots(f, X::IntervalLike{T}, ::Type{Bisection},
               strategy::SearchStrategy, tol::Float64) where {T}

    branch_and_prune(X, Bisection(f), strategy, tol)
end


# For `NewtonLike` acting on `Interval`
function _roots(f, X::Interval{T}, contractor::NewtonLike,
               strategy::SearchStrategy, tol::Float64) where {T}

    deriv = x -> ForwardDiff.derivative(f, x)
    _roots(f, deriv, X, contractor, strategy, tol)
end

function _roots(f, deriv, X::Interval{T}, contractor::NewtonLike,
               strategy::SearchStrategy, tol::Float64) where {T}

    branch_and_prune(X, contractor(f, deriv), strategy, tol)
end


# For `NewtonLike` acting on `IntervalBox`
function _roots(f, X::IntervalBox{T}, contractor::NewtonLike,
               strategy::SearchStrategy, tol::Float64) where {T}

    deriv = x -> ForwardDiff.jacobian(f, x)
    _roots(f, deriv, X, contractor, strategy, tol)
end

function _roots(f, deriv, X::IntervalBox{T}, contractor::NewtonLike,
               strategy::SearchStrategy, tol::Float64) where {T}

    branch_and_prune(X, contractor(f, deriv), strategy, tol)
end


# Acting on`Root`
# TODO: Use previous status information about roots:
function _roots(f, r::Root, contractor::Type{C},
               strategy::SearchStrategy, tol::Float64) where {C<:Contractor}

    _roots(f, r.interval, contractor, strategy, tol)
end

function _roots(f, deriv, r::Root, contractor::Type{C},
               strategy::SearchStrategy, tol::Float64) where {C<:Contractor}

    _roots(f, deriv, r.interval, contractor, strategy, tol)
end


# Acting on `Vector` of `Root`
function _roots(f, V::Vector{Root{T}}, contractor::Type{C},
               strategy::SearchStrategy, tol::Float64) where {T, C<:Contractor}

    vcat(_roots.(f, V, contractor, strategy, tol)...)
end

function _roots(f, deriv, V::Vector{Root{T}}, contractor::Type{C},
               strategy::SearchStrategy, tol::Float64) where {T, C<:Contractor}

    vcat(_roots.(f, deriv, V, contractor, strategy, tol)...)
end


# Acting on complex `Interval`
function _roots(f, Xc::Complex{Interval{T}}, contractor::Type{C},
               strategy::SearchStrategy, tol::Float64) where {T, C<:Contractor}

    g = realify(f)
    Y = IntervalBox(reim(Xc)...)
    rts = _roots(g, Y, contractor, strategy, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

function _roots(f, Xc::Complex{Interval{T}}, contractor::NewtonLike,
               strategy::SearchStrategy, tol::Float64) where {T}

    g = realify(f)
    g_prime = x -> ForwardDiff.jacobian(g, x)
    Y = IntervalBox(reim(Xc)...)
    rts = _roots(g, g_prime, Y, contractor, strategy, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

function _roots(f, deriv, Xc::Complex{Interval{T}}, contractor::NewtonLike,
               strategy::SearchStrategy, tol::Float64) where {T}

    g = realify(f)
    g_prime = realify_derivative(deriv)
    Y = IntervalBox(reim(Xc)...)
    rts = _roots(g, g_prime, Y, contractor, strategy, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end
