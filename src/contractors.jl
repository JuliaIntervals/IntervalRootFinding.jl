function image_contains_zero(f, R::Root)
    X = root_region(R)
    R.status == :empty && return Root(X, :empty)

    imX = f(X)

    if !(all(in_interval.(0, imX))) || isempty_region(imX)
        return Root(X, :empty)
    end

    return Root(X, :unknown)
end

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
    Y = inv(derivative(m))

    return m - Y*f(m) + (interval(1) - Y*derivative(X)) * (X - m)
end

function contract(::Type{Krawczyk}, f, derivative, X::AbstractVector)
    # We check the interval derivative first because if it works, we know that
    # the derivative of the center below works as well
    J = derivative(X)
    any(xx -> decoration(xx) == trv, J) && return interval(-Inf, Inf, trv)

    m = mid.(X)
    dm = derivative(m)
    det(dm) == 0 && return interval(-Inf, Inf, trv)

    Y = interval.(inv(dm))
    mm = interval.(m)

    return mm - Y*f(mm) + (interval(I) - Y*J) * (X - mm)
end

function contract(::Type{C}, f, derivative, R::Root) where {C <: AbstractContractor}
    # We first check with the simple bisection method
    # If we can prove it is empty at this point, we don't go further
    R2 = image_contains_zero(f, R)
    C == Bisection && return R2
    R2.status == :empty && return R2

    X = root_region(R)
    contracted_X = contract(C, f, derivative, X)
    istrivial(contracted_X) && return Root(X, :unknown)

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

function contract(root_problem::RootProblem{C}, R::Root) where C
    return contract(C, root_problem.f, root_problem.derivative, R)
end

"""
    refine(root_problem::RootProblem, X::Root)

Refine a root.
"""
function refine(root_problem::RootProblem{C}, R::Root) where C
    root_status(R) != :unique && throw(ArgumentError("trying to refine a non unique root"))

    X = root_region(R)

    while diam_region(X) > root_problem.abstol
        contracted = contract(C, root_problem.f, root_problem.derivative, X)
        istrivial(contracted) && break  # no point going further
        NX = intersect_region(contracted, X)
        isequal_region(NX, X) && break  # reached limit of precision
        X = NX
    end

    return Root(X, :unique)
end