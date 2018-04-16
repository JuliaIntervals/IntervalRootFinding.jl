
import IntervalArithmetic: diam, isinterior

export branch_and_prune, Bisection, Newton

diam(x::Root) = diam(x.interval)

Base.size(x::Interval) = (1,)

isinterior{N}(X::IntervalBox{N}, Y::IntervalBox{N}) = all(isinterior.(X, Y))

# contractors:
"""
    Contractor{F}

    Abstract type for contractors.
"""
abstract type Contractor{F} end

export Bisection, Newton

struct Bisection{F} <: Contractor{F}
    f::F
end

function (contractor::Bisection)(X)
    image = contractor.f(X)

    if !(contains_zero(image))
        return :empty, X
    end

    return :unknown, X
end


struct Newton{F,FP,O} <: Contractor{F}
    f::F
    f_prime::FP
    op::O
end

function Newton(f::Function, f_prime::Function)
    Newton(f, f_prime, N)
end

function Newton(f_prime::Function)
    return NewtonConstructor(f_prime)
end

function (C::Newton)(X)
    # use Bisection contractor for this:
    if !(contains_zero(IntervalBox(C.f(X))))
        return :empty, X
    end

    # given that have the Jacobian, can also do mean value form

    NX = C.op(C.f, C.f_prime, X) ∩ X

    isempty(NX) && return :empty, X

    if NX ⪽ X  # isinterior; know there's a unique root inside
        NX =  refine(X -> C.op(C.f, C.f_prime, X), NX)
        return :unique, NX
    end

    return :unknown, NX
end


struct NewtonConstructor{FP}
    f_prime::FP
end


"""
    branch_and_prune(X, contract, tol=1e-3)

Generic branch and prune routine for finding isolated roots using the `contract`
function as the contractor.

Inputs:
- `X`: `Interval` or `IntervalBox`
- `contractor`: function that determines the status of a given box `X`. It
    returns the new box and a symbol indicating the status. Current possible
    values are of type `Bisection` or `Newton`

"""
function branch_and_prune(X, contractor, tol=1e-3)
    working = [X]
    outputs = Root{typeof(X)}[]

    sizehint!(outputs, 100)
    sizehint!(working, 1000)

    while !isempty(working)
        # @show working
        X = pop!(working)

        status, output = contractor(X)

        if status == :empty
            continue

        elseif status == :unique
            push!(outputs, Root(output, :unique))

        elseif diam(output) < tol
            push!(outputs, Root(output, :unknown))

        else  # branch
            X1, X2 = bisect(X)

            push!(working, X1, X2)
        end
    end

    return outputs
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


contains_zero(X::Interval{T}) where {T} = zero(T) ∈ X
contains_zero(X::SVector) = all(contains_zero.(X))
contains_zero(X::IntervalBox) = all(contains_zero.(X))


"""
Generic refine operation for Krawczyk and Newton.
This assumes that it is already known that `X` contains a unique root.
Call using e.g. `op = X -> N(f, f_prime, X)`
"""
function refine(op, X)

    tolerance = 1e-16

    while diam(X) > tolerance  # avoid problem with tiny floating-point numbers if 0 is a root
        NX = op(X) ∩ X
        NX == X && break  # reached limit of precision
        X = NX
    end

    return X
end


IntervalLike{T} = Union{Interval{T}, IntervalBox{T}}

"""
    roots(f, X, contractor, tol=1e-3)

Generic branch and prune routine for finding isolated roots of a function
`f:R^n → R^n` in a box `X`, or a vector of boxes.

Inputs:
- `f`: function whose roots will be found
- `X`: `Interval` or `IntervalBox`
- `contractor`: function that, when applied to the function `f`, determines
    the status of a given box `X`. It returns the new box and a symbol indicating
    the status. Current possible values are `Bisection`, `Newton` and
    `Newton(f_prime)` where `f_prime` is the derivative or jacobian of `f`.

"""
# Contractor specific `roots` functions
function roots(f, X::IntervalLike{T}, ::Type{Bisection}, tol::Float64=1e-3) where {T}
    branch_and_prune(X, Bisection(f), tol)
end

function roots(f, X::Interval{T}, ::Type{Newton}, tol::Float64=1e-3) where {T}
    branch_and_prune(X, Newton(f, x -> ForwardDiff.derivative(f, x)), tol)
end

function roots(f, X::IntervalBox{T}, ::Type{Newton}, tol::Float64=1e-3) where {T}
    branch_and_prune(X, Newton(f, x -> ForwardDiff.jacobian(f, x)), tol)
end

function roots{T}(f, X::IntervalLike{T}, nc::NewtonConstructor, tol::Float64=1e-3)
    branch_and_prune(X, Newton(f, nc.f_prime), tol)
end

# `roots` function for cases where `X` is not an `Interval` or `IntervalBox`
function roots(f, V::Vector{Root{T}}, contractor::Type{C}, tol::Float64=1e-3) where {T, C<:Contractor}
    reduce(append!, Root{T}[], [roots(f, X.interval, contractor, tol) for X in V])
end

function roots(f, V::Vector{T}, contractor::Type{C}, tol::Float64=1e-3) where {T, C<:Contractor}
    reduce(append!, Root{T}[], [roots(f, X, contractor, tol) for X in V])
end

function roots(f, Xc::Complex{Interval{T}}, contractor::Type{C}, tol::Float64=1e-3) where {T, C<:Contractor}
    g = realify(f)
    Y = IntervalBox(reim(Xc))
    rts = roots(g, Y, contractor, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

function roots(f, Xc::Complex{Interval{T}}, nc::NewtonConstructor, tol::Float64=1e-3) where {T}
    g = realify(f)
    g_prime = realify_derivative(nc.f_prime)
    Y = IntervalBox(reim(Xc))
    rts = roots(g, Y, Newton(g_prime), tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

# Default
roots(f, X, tol::Float64=1e-3) = roots(f, X, Newton, tol)
