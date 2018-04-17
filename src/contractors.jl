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
    image = contractor.f(X)

    if !(contains_zero(image))
        return :empty, X
    end

    return :unknown, X
end

# Newton
struct Newton{F,FP,O} <: Contractor{F}
    f::F
    f_prime::FP
    op::O
end

function Newton(f::Function, f_prime::Function)
    Newton(f, f_prime, 𝒩)
end

function (C::Newton)(X, tol)
    # use Bisection contractor for this:
    if !(contains_zero(IntervalBox(C.f(X))))
        return :empty, X
    end

    # given that have the Jacobian, can also do mean value form

    NX = C.op(C.f, C.f_prime, X) ∩ X

    isempty(NX) && return :empty, X

    if NX ⪽ X  # isinterior; know there's a unique root inside
        NX =  refine(X -> C.op(C.f, C.f_prime, X), NX, tol)
        return :unique, NX
    end

    return :unknown, NX
end


# What is this guarded_mid for? Shouldn't it be checking if f(m)==0?
doc"""Returns the midpoint of the interval x, slightly shifted in case
the midpoint is an exact root"""
function guarded_mid{T}(f, x::Interval{T})
    m = mid(x)

    if f(m) == 0                      # midpoint exactly a root
        α = convert(T, 0.46875)      # close to 0.5, but exactly representable as a floating point
        m = α*x.lo + (one(T)-α)*x.hi   # displace to another point in the interval
    end

    @assert m ∈ x

    m
end

doc"""
Single-variable Newton operator
"""
function 𝒩{T}(f, x::Interval{T}, deriv::Interval{T})
    m = Interval( guarded_mid(f, x) )

    m - (f(m) / deriv)
end

function 𝒩{T}(f, x::Interval{T})
    m = Interval( guarded_mid(f, x) )

    m - (f(m) / ForwardDiff.derivative(f, x))
end

function 𝒩{T}(f, f_prime, X::Interval{T})
    m = Interval( guarded_mid(f, X) )

    m - (f(m) / f_prime(X))
end



doc"""
Multi-variable Newton operator.
Requires the function to be defined using the `@intervalbox` macro.
"""
function 𝒩(f::Function, jacobian::Function, X::IntervalBox)  # multidimensional Newton operator

    m = IntervalBox(Interval.(mid(X)))
    J = jacobian(SVector(X))

    return IntervalBox( (m - (J \ f(m))... ) )
end



"""
Generic refine operation for Krawczyk and Newton.
This function assumes that it is already known that `X` contains a unique root.
Call using e.g. `op = X -> N(f, f_prime, X)`
"""
function refine(op, X, tolerance=1e-16)

    while diam(X) > tolerance  # avoid problem with tiny floating-point numbers if 0 is a root
        NX = op(X) ∩ X
        NX == X && break  # reached limit of precision
        X = NX
    end

    return X
end
