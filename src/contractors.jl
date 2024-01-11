export Bisection, Newton, Krawczyk

"""
    AbstractContractor

Abstract type for contractors.
"""
abstract type AbstractContractor end

struct Bisection <: AbstractContractor end
struct Newton <: AbstractContractor end
struct Krawczyk <: AbstractContractor end

function contract(::Type{Newton}, f, derivative, X::Interval, where_mid)
    m = interval(mid(X, where_mid))
    return m - (f(m) / derivative(X))
end

function contract(::Type{Newton}, f, derivative, X::AbstractVector, where_mid)
    m = interval.(mid.(X, where_mid))
    J = derivative(X)
    y = gauss_elimination_interval(J, f(m))  # J \ f(m)
    return m .- y
end

function contract(::Type{Krawczyk}, f, derivative, X::Interval, where_mid)
    m = interval(mid(X, where_mid))
    Y = 1 / derivative(m)

    return m - Y*f(m) + (1 - Y*derivative(X)) * (X - m)
end

function contract(::Type{Krawczyk}, f, derivative, X::AbstractVector, where_mid)
    mm = mid.(X, where_mid)
    J = derivative(X)
    Y = mid.(inv(derivative(mm)))

    return m - Y*f(mm) + (I - Y*J) * (X.v - m)
end

function contract(root_problem::RootProblem{C}, X::region) where C
    CX = contract(C, root_problem.f, root_problem.derivative, region(X), root_problem.where_bisect)
    return Region(CX)
end

"""
    safe_isempty(X)

Similar to `isempty` function for `IntervalBox`, but also works for `SVector`
of `Interval`.
"""
safe_isempty(X) = any(isempty_interval.(X))

function image_contains_zero(f, R::Root)
    X = interval(R)
    R.status == :empty && return Root(X, :empty)

    imX = f(X)

    if !(all(in_interval.(0, imX))) || safe_isempty(imX)
        return Root(X, :empty)
    end

    return Root(X, :unknown)
end

function contract_root(root_problem::RootProblem{C}, R::Root) where C
    # We first check with the simple bisection method
    # If we can prove it is empty at this point, we don't go further
    R2 = image_contains_zero(root_problem.f, R)
    C == Bisection && return R2
    R2.status == :empty && return R2

    X = interval(R)
    contracted_X = contract(root_problem, X)

    # Only happens if X is partially out of the domain of f
    safe_isempty(contracted_X) && return Root(X, :unknown)  # force bisection

    NX = intersect_interval(bareinterval(contracted_X), bareinterval(X))
    NX = IntervalArithmetic._unsafe_interval(NX, min(decoration(contracted_X), decoration(X)), isguaranteed(contracted_X))

    !isbounded(X) && return Root(NX, :unknown)  # force bisection
    safe_isempty(NX) && return Root(X, :empty)

    if R.status == :unique || NX ⪽ X  # isstrictsubset_interval, we know there's a unique root inside
        return Root(NX, :unique)
    end

    return Root(NX, :unknown)
end

"""
    refine(op, X::Root, tol)

Wrap the refine method to leave unchanged intervals that are not guaranteed to
contain an unique solution.
"""
function refine(root_problem::RootProblem, R::Root)
    root_status(R) != :unique && return R
    return Root(refine_root(root_problem, region(R)))
end

"""
    refine(C, X::Region, tol)

Refine a interval known to contain a solution.

This function assumes that it is already known that `X` contains a unique root.
"""
function refine_root(root_problem::RootProblem, X::Region)
    while diam(X) > root_problem.abstol
        NX = intersect_interval(bareinterval(C(X)), bareinterval(X))
        NX = IntervalArithmetic._unsafe_interval(NX, min(decoration(C(X)), decoration(X)), isguaranteed(C(X)))
        isequal_interval(NX, X) && break  # reached limit of precision
        X = NX
    end

    return X
end