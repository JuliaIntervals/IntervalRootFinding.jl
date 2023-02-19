export Contractor
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

"""
    Newton{F, FP} <: Contractor{F}

Contractor type for the interval Newton method.

# Fields
    - `f::F`: function whose roots are searched
    - `f::FP`: derivative or jacobian of `f`

-----

    (N::Newton)(X, α=where_bisect)

Contract an interval `X` using Newton operator and return the
contracted interval together with its status.

# Inputs
    - `R`: Root object containing the interval to contract.
    - `α`: Point of bisection of intervals.
"""
struct Newton{F, FP} <: Contractor{F}
    f::F
    f′::FP   # use \prime<TAB> for ′
end

function (N::Newton)(X::Interval ; α=where_bisect)
    m = Interval(mid(X, α))
    return m - (N.f(m) / N.f′(X))
end

function (N::Newton)(X::IntervalBox ; α=where_bisect)
    m = Interval.(mid(X, α))
    J = N.f′(X)
    y = gauss_elimination_interval(J, N.f(m))  # J \ f(m)
    return IntervalBox(m .- y)
end

"""
    Krawczyk{F, FP} <: Contractor{F}

Contractor type for the interval Krawczyk method.

# Fields
    - `f::F`: function whose roots are searched
    - `f::FP`: derivative or jacobian of `f`

-----

    (K::Krawczyk)(X ; α=where_bisect)

Contract an interval `X` using Krawczyk operator and return the
contracted interval together with its status.

# Inputs
    - `R`: Root object containing the interval to contract.
    - `α`: Point of bisection of intervals.
"""
struct Krawczyk{F, FP} <: Contractor{F}
    f::F
    f′::FP   # use \prime<TAB> for ′
end

function (K::Krawczyk)(X::Interval ; α=where_bisect)
    m = Interval(mid(X, α))
    Y = 1 / K.f′(m)

    return m - Y*K.f(m) + (1 - Y*K.f′(X)) * (X - m)
end

function (K::Krawczyk)(X::IntervalBox ; α=where_bisect)
    jacobian = K.f′
    m = mid(X, α)
    mm = IntervalBox(m)
    J = jacobian(X)
    Y = mid.(inv(jacobian(mm)))

    return m - Y*K.f(mm) + (I - Y*J) * (X.v - m)
end

"""
    safe_isempty(X)

Similar to `isempty` function for `IntervalBox`, but also works for `SVector`
of `Interval`.
"""
safe_isempty(X) = isempty(IntervalBox(X))

"""
    contract(contractor, R)

Contraction operation for contractors using the first derivative of the
function.
"""
function contract(B::Bisection, R::Root)
    X = interval(R)
    R.status == :empty && return Root(X, :empty)

    imX = B.f(X)

    if !(contains_zero(imX)) || safe_isempty(imX)
        return Root(X, :empty)
    end

    return Root(X, :unknown)
end

function contract(C::Union{Newton, Krawczyk}, R::Root)
    # We first check with the simple bisection method
    # If we can prove it is empty at this point, we don't go further
    R2 = contract(Bisection(C.f), R)
    R2.status == :empty && return R2

    X = interval(R)
    contracted_X = C(X)

    # Only happens if X is partially out of the domain of f
    safe_isempty(contracted_X) && return Root(X, :unknown)  # force bisection

    NX = contracted_X ∩ X

    isinf(X) && return Root(NX, :unknown)  # force bisection
    safe_isempty(NX) && return Root(X, :empty)

    if R.status == :unique || NX ⪽ X  # isinterior, we know there's a unique root inside
        return Root(NX, :unique)
    end

    return Root(NX, :unknown)
end

"""
    refine(C, X::Region, tol)

Refine a interval known to contain a solution.

This function assumes that it is already known that `X` contains a unique root.
"""
function refine(C::Union{Newton, Krawczyk}, X::Region, root_problem)
    while diam(X) > root_problem.abstol
        NX = C(X) ∩ X
        NX == X && break  # reached limit of precision
        X = NX
    end

    return X
end

"""
    refine(op, X::Tuple{Symbol, Region}, tol)

Wrap the refine method to leave unchanged intervals that are not guaranteed to
contain an unique solution.
"""
function refine(C::Contractor, R::Root, root_problem)
    root_status(R) != :unique && return R
    return Root(refine(C, interval(R), root_problem), :unique)
end
