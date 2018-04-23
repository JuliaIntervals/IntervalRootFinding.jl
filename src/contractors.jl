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

Base.isinf(X::IntervalBox) = any(isinf.(X))

function (C::Newton)(X, tol)
    # use Bisection contractor for this:
    if !(contains_zero(C.f(X)))
        return :empty, X
    end

    # given that have the Jacobian, can also do mean value form

    NX = 𝒩(C.f, C.f′, X) ∩ X

    isempty(NX) && return :empty, X

    if isinf(X)
        return :unknown, NX  # force bisection
    end

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


IntervalArithmetic.mid(X::IntervalBox, α) = mid.(X, α)

doc"""
Multi-variable Newton operator.
"""
function 𝒩(f::Function, jacobian::Function, X::IntervalBox)  # multidimensional Newton operator

    m = IntervalBox(Interval.(mid(X, where_bisect)))
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
