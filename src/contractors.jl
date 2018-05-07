export Bisection, Newton, Krawczyk

Base.isinf(X::IntervalBox) = any(isinf.(X))
IntervalArithmetic.mid(X::IntervalBox, α) = mid.(X, α)

doc"""
    Contractor{F}

    Abstract type for contractors.
"""
abstract type Contractor{F} end

doc"""
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

doc"""
    newtonlike_contract(op, X, tol)

    Contraction operation for contractors using the first derivative of the
    function. This contraction use a bisection scheme to refine the intervals
    with `:unkown` status.

    Currently `Newton` and `Krawczyk` contractors uses this.
"""
function newtonlike_contract(op, C, X, tol)
    # use Bisection contractor for this:
    if !(contains_zero(C.f(X)))
        return :empty, X
    end

    # given that have the Jacobian, can also do mean value form

    NX = op(C.f, C.f′, X) ∩ X

    isempty(NX) && return :empty, X
    isinf(X) && return :unknown, NX  # force bisection

    if NX ⪽ X  # isinterior; know there's a unique root inside
        NX =  refine(X -> op(C.f, C.f′, X), NX, tol)
        return :unique, NX
    end

    return :unknown, NX
end

doc"""
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


doc"""
Single-variable Newton operator
"""
function 𝒩{T}(f, X::Interval{T})
    m = Interval(mid(X, where_bisect))

    m - (f(m) / ForwardDiff.derivative(f, X))
end

function 𝒩{T}(f, f′, X::Interval{T})
    m = Interval(mid(X, where_bisect))

    m - (f(m) / f′(X))
end

function 𝒩{T}(f, X::Interval{T}, dX::Interval{T})
    m = Interval(mid(X, where_bisect))

    m - (f(m) / dX)
end

doc"""
Multi-variable Newton operator.
"""
function 𝒩(f::Function, jacobian::Function, X::IntervalBox)  # multidimensional Newton operator
    m = IntervalBox(Interval.(mid(X, where_bisect)))
    J = jacobian(SVector(X))

    return IntervalBox(m - (J \ f(m)))
end


doc"""
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


doc"""
Single-variable Krawczyk operator
"""
function 𝒦(f, f′, X::Interval{T}) where {T}
    m = Interval(mid(X))
    Y = 1/f′(m)

    m - Y*f(m) + (1 - Y*f′(X))*(X - m)
end

doc"""
Multi-variable Krawczyk operator
"""
function 𝒦(f, jacobian, X::IntervalBox{T}) where {T}
    m = mid(X)
    J = jacobian(X)
    Y = inv(jacobian(m))
    m = IntervalBox(Interval.(m))

    IntervalBox(m - Y*f(m) + (I - Y*J)*(X - m))
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
