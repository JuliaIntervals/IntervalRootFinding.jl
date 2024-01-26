export Bisection, Newton, Krawczyk

"""
    AbstractContractor

Abstract type for contractors.
"""
abstract type AbstractContractor end

struct Bisection <: AbstractContractor end
struct Newton <: AbstractContractor end
struct Krawczyk <: AbstractContractor end

function contract(::Type{Newton}, f, derivative, X::Interval)
    m = interval(mid(X))
    return m - (f(m) / derivative(X))
end

function contract(::Type{Newton}, f, derivative, X::AbstractVector)
    m = interval.(mid.(X))
    J = derivative(X)
    y = gauss_elimination_interval(J, f(m))  # J \ f(m)
    return m .- y
end

function contract(::Type{Krawczyk}, f, derivative, X::Interval)
    m = interval(mid(X))
    Y = 1 / derivative(m)

    return m - Y*f(m) + (1 - Y*derivative(X)) * (X - m)
end

function contract(::Type{Krawczyk}, f, derivative, X::AbstractVector)
    m = mid.(X)

    dm = derivative(m)
    det(dm) == 0 && return interval(-Inf, Inf, trv)

    Y = inv(dm)
    J = derivative(X)
    mm = interval.(m)

    return mm - Y*f(mm) + (I - Y*J) * (X - mm)
end

function contract(root_problem::RootProblem{C}, X) where C
    contracted = contract(C, root_problem.f, root_problem.derivative, X)
    istrivial(contracted) && return X
    return contracted
end

function image_contains_zero(f, R::Root)
    X = root_region(R)
    R.status == :empty && return Root(X, :empty)

    imX = f(X)

    if !(all(in_interval.(0, imX))) || isempty_region(imX)
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

    X = root_region(R)
    contracted_X = contract(root_problem, X)

    # Only happens if X is partially out of the domain of f
    isempty_region(contracted_X) && return Root(X, :unknown)  # force bisection

    NX = intersect_region(contracted_X, X)

    !isbounded_region(X) && return Root(NX, :unknown)  # force bisection
    isempty_region(NX) && return Root(X, :empty)

    if R.status == :unique || NX ⪽ X  # isstrictsubset_interval, we know there's a unique root inside
        return Root(NX, :unique)
    end

    return Root(NX, :unknown)
end

"""
    refine(root_problem::RootProblem, X::Root)

Refine a root.
"""
function refine(root_problem::RootProblem, R::Root)
    root_status(R) != :unique && return R

    X = root_region(R)

    while diam_region(X) > root_problem.abstol
        NX = intersect_region(contract(root_problem, X), X)
        isequal_region(NX, X) && break  # reached limit of precision
        X = NX
    end

    return Root(X, :unique)
end