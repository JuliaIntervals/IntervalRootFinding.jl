
import IntervalArithmetic: diam, isinterior

export branch_and_prune, Bisection, Newton

diam(x::Root) = diam(x.interval)

"""
    branch_and_prune(X, contractor, tol=1e-3)

Generic branch and prune routine for finding isolated roots using the `contract`
function as the contractor.

Inputs:
- `X`: `Interval` or `IntervalBox`
- `contractor`: function that determines the status of a given box `X`. It
    returns the new box and a symbol indicating the status. Current possible
    values are of type `Bisection` or `Newton`

"""
function branch_and_prune(X, contractor, tol=1e-3)
    iter = RootSearch(X, contractor, tol)
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
    roots(f, X, contractor, tol=1e-3)
    roots(f, deriv, X, contractor, tol=1e-3)

Uses a generic branch and prune routine to find in principle all isolated roots of a function
`f:R^n → R^n` in a box `X`, or a vector of boxes.

Inputs:
- `f`: function whose roots will be found
- `X`: `Interval` or `IntervalBox`
- `contractor`: function that, when applied to the function `f`, determines
    the status of a given box `X`. It returns the new box and a symbol indicating
    the status. Current possible values are `Bisection`, `Newton` and `Krawczyk`
- `deriv`: explicit derivative of `f` for `Newton` and `Krawczyk`

"""
# Contractor specific `roots` functions
function roots(f, X::IntervalLike{T}, ::Type{Bisection}, tol::Float64=1e-3) where {T}
    branch_and_prune(X, Bisection(f), tol)
end

function roots(f, X::Interval{T}, C::NewtonLike, tol::Float64=1e-3) where {T}

    deriv = x -> ForwardDiff.derivative(f, x)

    roots(f, deriv, X, C, tol)
end

function roots(f, deriv, X::Interval{T}, C::NewtonLike, tol::Float64=1e-3) where {T}
    branch_and_prune(X, C(f, deriv), tol)
end

function roots(f, X::IntervalBox{T}, C::NewtonLike, tol::Float64=1e-3) where {T}

    deriv = x -> ForwardDiff.jacobian(f, x)

    roots(f, deriv, X, C, tol)
end

function roots(f, deriv, X::IntervalBox{T}, C::NewtonLike, tol::Float64=1e-3) where {T}
    branch_and_prune(X, C(f, deriv), tol)
end

roots(f, r::Root, contractor::Type{C}, tol::Float64=1e-3) where {C<:Contractor} =
    roots(f, r.interval, contractor, tol)
roots(f, deriv, r::Root, contractor::Type{C}, tol::Float64=1e-3) where {C<:Contractor} =
    roots(f, deriv, r.interval, contractor, tol)


# Acting on a Vector:

# TODO: Use previous status information about roots:
roots(f, V::Vector{Root{T}}, contractor::Type{C}, tol::Float64=1e-3) where {T, C<:Contractor} =
    vcat(roots.(f, V, contractor, tol)...)
roots(f, deriv, V::Vector{Root{T}}, contractor::Type{C}, tol::Float64=1e-3) where {T, C<:Contractor} =
    vcat(roots.(f, deriv, V, contractor, tol)...)



# Complex:

function roots(f, Xc::Complex{Interval{T}}, contractor::Type{C},
        tol::Float64=1e-3) where {T, C<:Contractor}

    g = realify(f)
    Y = IntervalBox(reim(Xc)...)
    rts = roots(g, Y, contractor, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

function roots(f, Xc::Complex{Interval{T}}, C::NewtonLike, tol::Float64=1e-3) where {T}

    g = realify(f)

    g_prime = x -> ForwardDiff.jacobian(g, x)

    Y = IntervalBox(reim(Xc)...)
    rts = roots(g, g_prime, Y, C, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

function roots(f, deriv, Xc::Complex{Interval{T}}, C::NewtonLike, tol::Float64=1e-3) where {T}

    g = realify(f)

    g_prime = realify_derivative(deriv)

    Y = IntervalBox(reim(Xc)...)
    rts = roots(g, g_prime, Y, C, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

# Default
roots(f::F, deriv::F, X, tol::Float64=1e-15) where {F<:Function} = roots(f, deriv, X, Newton, tol)
roots(f::F, X, tol::Float64=1e-15) where {F<:Function} = roots(f, X, Newton, tol)
