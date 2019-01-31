export Bisection, Newton, Krawczyk

"""
    𝒩(f, f′, X, α)

Single-variable Newton operator.

The symbol for the operator is accessed with `\\scrN<tab>`.
"""
function 𝒩(f, f′, X::Interval{T}, α) where {T}
    m = Interval(mid(X, α))

    return m - (f(m) / f′(X))
end

"""
    𝒩(f, jacobian, X, α)

Multi-variable Newton operator.
"""
function 𝒩(f::Function, jacobian::Function, X::IntervalBox, α)  # multidimensional Newton operator
    m = IntervalBox(Interval.(mid(X, α)))
    J = jacobian(X)

    return IntervalBox(m .- (J \ f(m)))
end

"""
    𝒦(f, f′, X, α)

Single-variable Krawczyk operator.

The symbol for the operator is accessed with `\\scrK<tab>`.
"""
function 𝒦(f, f′, X::Interval{T}, α) where {T}
    m = Interval(mid(X, α))
    Y = 1 / f′(m)

    return m - Y*f(m) + (1 - Y*f′(X)) * (X - m)
end

"""
    𝒦(f, jacobian, X, α)

Multi-variable Krawczyk operator.
"""
function 𝒦(f, jacobian, X::IntervalBox{T}, α) where {T}
    m = mid(X, α)
    mm = IntervalBox(m)
    J = jacobian(X)
    Y = mid.(inv(jacobian(mm)))

    return m - Y*f(mm) + (I - Y*J) * (X.v - m)
end


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

for (Method, 𝒪) in ((:Newton, 𝒩), (:Krawczyk, 𝒦))
    doc = """
        $Method{F, FP} <: Contractor{F}

    Contractor type for the $Method method.

    # Fields
        - `f::F`: function whose roots are searched
        - `f::FP`: derivative or jacobian of `f`

    -----

        (C::$Method)(X, tol, α=where_bisect)

    Contract an interval `X` using $Method operator and return the
    contracted interval together with its status.

    # Inputs
        - `X`: Interval to contract.
        - `tol`: Precision to which unique solutions are refined.
        - `α`: Point of bisection of intervals.
    """

    @eval begin
        struct $Method{F, FP} <: Contractor{F}
            f::F
            f′::FP   # use \prime<TAB> for ′
        end

        function (C::$Method)(X, tol, α=where_bisect)
            op = x -> $𝒪(C.f, C.f′, x, α)
            rt = determine_region_status(op, C.f, X)
            return refine(op, rt, tol)
        end

        @doc $doc $Method
    end
end


"""
    safe_isempty(X)

Similar to `isempty` function for `IntervalBox`, but also works for `SVector`
of `Interval`.
"""
safe_isempty(X) = isempty(IntervalBox(X))

"""
    determine_region_status(contract, f, X)

Contraction operation for contractors using the first derivative of the
function.

Currently `Newton` and `Krawczyk` contractors use this.
"""
function determine_region_status(op, f, X)
    imX = f(X)

    !(contains_zero(imX)) && return :empty, X

    safe_isempty(imX) && return :empty, X  # X is fully outside of the domain of f

    contracted_X = op(X)

    # Only happens if X is partially out of the domain of f
    safe_isempty(contracted_X) && return :unknown, X  # force bisection

    # given that have the Jacobian, can also do mean value form
    NX = contracted_X ∩ X

    isinf(X) && return :unknown, NX  # force bisection
    safe_isempty(NX) && return :empty, X

    NX ⪽ X  && return :unique, NX  # isinterior; know there's a unique root inside

    return :unknown, NX
end

"""
    refine(op, X::Region, tol)

Generic refine operation for Krawczyk and Newton.
This function assumes that it is already known that `X` contains a unique root.
"""
function refine(op, X::Region, tol)
    while diam(X) > tol  # avoid problem with tiny floating-point numbers if 0 is a root
        NX = op(X) ∩ X
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
function refine(op, rt::Tuple{Symbol, Region}, tol)
    rt[1] != :unique && return rt
    return :unique, refine(op, rt[2], tol)
end
