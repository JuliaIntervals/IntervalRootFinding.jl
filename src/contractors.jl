# contractors:
"""
    Contractor{F}

    Abstract type for contractors.
"""
abstract type Contractor{F} end

export Bisection, Newton

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

function (C::Newton)(X, tol)
    # use Bisection contractor for this:
    if !(contains_zero(C.f(X)))
        return :empty, X
    end

    # given that have the Jacobian, can also do mean value form

    NX = 𝒩(C.f, C.f′, X) ∩ X

    isempty(NX) && return :empty, X

    if NX ⪽ X  # isinterior; know there's a unique root inside
        NX =  refine(X -> 𝒩(C.f, C.f′, X), NX, tol)
        return :unique, NX
    end

    return :unknown, NX
end



doc"""
Single-variable Newton operator
"""
function 𝒩{T}(f, x::Interval{T})
    m = Interval(mid(X))

    m - (f(m) / ForwardDiff.derivative(f, x))
end

function 𝒩{T}(f, f′, X::Interval{T})
    m = Interval(mid(X))

    m - (f(m) / f′(X))
end



doc"""
Multi-variable Newton operator.
"""
function 𝒩(f::Function, jacobian::Function, X::IntervalBox)  # multidimensional Newton operator

    m = IntervalBox(Interval.(mid(X)))
    J = jacobian(SVector(X))

    return IntervalBox(m - (J \ f(m)))
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
