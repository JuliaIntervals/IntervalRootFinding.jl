export Bisection, Newton, Krawczyk


"""
    Contractor{F}

    Abstract type for contractors.
"""
abstract type Contractor{F} end

"""
    Bisection{F} <: Contractor{F}

    Contractor type for the bisection method.
"""
struct Bisection{F} <: Contractor{F}
    f::F
end

function (contractor::Bisection)(X, tol)
    image = (contractor.f)(X)

    if !(contains_zero(image))
        return :empty, X
    end

    return :unknown, X
end

"""
    newtonlike_contract(op, X, tol)

    Contraction operation for contractors using the first derivative of the
    function. This contraction use a bisection scheme to refine the intervals
    with `:unkown` status.

    Currently `Newton` and `Krawczyk` contractors uses this.
"""
function newtonlike_contract(op, C, X, tol)
    imX = C.f(X)
    !(contains_zero(imX)) && return :empty, X
    # Only happens if X is fully outside the domain of f
    isempty(imX) && return :empty, X

    contracted_X = op(C.f, C.f′, X)

    # Only happens if X is partially out of the domain of f
    isempty(contracted_X) && return :unkown, X  # force bisection

    # given that have the Jacobian, can also do mean value form
    NX = contracted_X ∩ X

    isinf(X) && return :unknown, NX  # force bisection
    isempty(NX) && return :empty, X

    if NX ⪽ X  # isinterior; know there's a unique root inside
        NX =  refine(X -> op(C.f, C.f′, X), NX, tol)
        return :unique, NX
    end

    return :unknown, NX
end

"""
    Newton{F, FP} <: Contractor{F}

    Contractor type for the Newton method.

    # Fields
        - `f::F`: function whose roots are searched
        - `f::FP`: derivative or jacobian of `f`
"""
struct Newton{F,FP} <: Contractor{F}
    f::F
    f′::FP   # use \prime<TAB> for ′
end

function (C::Newton)(X, tol)
    newtonlike_contract(𝒩, C, X, tol)
end


"""
Single-variable Newton operator
"""
function 𝒩(f, X::Interval{T}) where {T}
    m = Interval(mid(X, where_bisect))

    m - (f(m) / ForwardDiff.derivative(f, X))
end

function 𝒩(f, f′, X::Interval{T}) where {T}
    m = Interval(mid(X, where_bisect))

    m - (f(m) / f′(X))
end

function 𝒩(f, X::Interval{T}, dX::Interval{T}) where {T}
    m = Interval(mid(X, where_bisect))

    m - (f(m) / dX)
end

"""
Multi-variable Newton operator.
"""
function 𝒩(f::Function, jacobian::Function, X::IntervalBox)  # multidimensional Newton operator
    m = IntervalBox(Interval.(mid(X, where_bisect)))
    J = jacobian(X)

    return IntervalBox(m .- (J \ f(m)))
end


"""
    Krawczyk{F, FP} <: Contractor{F}

    Contractor type for the Krawczyk method.

    # Fields
        - `f::F`: function whose roots are searched
        - `f::FP`: derivative or jacobian of `f`
"""
struct Krawczyk{F, FP} <: Contractor{F}
    f::F
    f′::FP   # use \prime<TAB> for ′
end

function (C::Krawczyk)(X, tol)
    newtonlike_contract(𝒦, C, X, tol)
end


"""
Single-variable Krawczyk operator
"""
function 𝒦(f, f′, X::Interval{T}) where {T}
    m = Interval(mid(X))
    Y = 1 / f′(m)

    m - Y*f(m) + (1 - Y*f′(X)) * (X - m)
end

"""
Multi-variable Krawczyk operator
"""
function 𝒦(f, jacobian, X::IntervalBox{T}) where {T}
    m = mid(X)
    J = jacobian(X)
    Y = inv(jacobian(m))
    mm = IntervalBox(m)

    res = m - Y*f(mm) + (I - Y*J) * (X.v - m)    # IntervalBox(res)
end

"""
Generic refine operation for Krawczyk and Newton.
This function assumes that it is already known that `X` contains a unique root.
Call using e.g. `op = X -> N(f, f′, X)`
"""
function refine(op, X, tolerance=1e-16)

    while diam(X) > tolerance  # avoid problem with tiny floating-point numbers if 0 is a root
        NX = op(X) ∩ X
        NX == X && break  # reached limit of precision
        X = NX
    end

    return X
end
