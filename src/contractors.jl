# contractors:
"""
    Contractor{F}

    Abstract type for contractors.
"""
abstract type Contractor{F} end

export Bisection, Newton, Krawczyk

# bisection
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

# Newton
struct Newton{F,FP} <: Contractor{F}
    f::F
    f′::FP   # use \prime<TAB> for ′
end

Base.isinf(X::IntervalBox) = any(isinf.(X))

function (C::Newton)(X, tol)
    # use Bisection contractor for this:
    if !(contains_zero(C.f(X)))
        return :empty, X
    end

    # given that have the Jacobian, can also do mean value form

    NX = 𝒩(C.f, C.f′, X) ∩ X

    isempty(NX) && return :empty, X
    isinf(X) && return :unknown, NX  # force bisection

    if NX ⪽ X  # isinterior; know there's a unique root inside
        NX =  refine(X -> 𝒩(C.f, C.f′, X), NX, tol)
        return :unique, NX
    end

    return :unknown, NX
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


doc"""
Multi-variable Newton operator.
"""
function 𝒩(f::Function, jacobian::Function, X::IntervalBox)  # multidimensional Newton operator
    m = IntervalBox(Interval.(mid(X, where_bisect)))
    J = jacobian(SVector(X))

    return IntervalBox(m - (J \ f(m)))
end


# Krawczyk
struct Krawczyk{F, FP} <: Contractor{F}
    f::F
    f′::FP   # use \prime<TAB> for ′
end

function (C::Krawczyk)(X, tol)
    # use Bisection contractor for this:
    if !(contains_zero(C.f(X)))
        return :empty, X
    end

    KX = 𝒦(C.f, C.f′, X) ∩ X

    isempty(KX) && return :empty, X
    isinf(X) && return :unknown, KX  # force bisection

    if KX ⪽ X  # isinterior; know there's a unique root inside
        KX =  refine(X -> 𝒦(C.f, C.f′, X), KX, tol)
        return :unique, KX
    end

    return :unknown, KX
end


doc"""
Single-variable Krawczyk operator
"""
function 𝒦(f, f′, X::Interval{T}) where {T}
    m = Interval(mid(X, where_bisect))
    Y = 1/f′(m)

    m - Y*f(m) + (1 - Y*f′(X))*(X - m)
end


@generated function 𝒦(f, jacobian, X::IntervalBox{N, T}) where {N, T}
    unit = eye(N)

    ex = quote
        m = mid(X, where_bisect)
        J = jacobian(X)
        Y = inv(jacobian(m))
        m = IntervalBox(Interval.(m))

        IntervalBox(m - Y*f(m) + ($unit - Y*J)*(X - m))
    end
    ex
end

IntervalArithmetic.mid(X::IntervalBox, α) = mid.(X, α)

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
