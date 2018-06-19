# This file is part of the ValidatedNumerics.jl package; MIT licensed

# Krawczyk method, following Tucker

# const D = differentiate

doc"""Returns two intervals, the first being a point within the
interval x such that the interval corresponding to the derivative of f there
does not contain zero, and the second is the inverse of its derivative"""
function guarded_derivative_midpoint(f::Function, f_prime::Function, x::Interval{T}) where {T}

    α = convert(T, 0.46875)   # close to 0.5, but exactly representable as a floating point

    m = Interval(mid(x))

    C = inv(f_prime(m))

    # Check that 0 is not in C; if so, consider another point rather than m
    i = 0
    while zero(T) ∈ C || isempty(C)
        m = Interval( α*x.lo + (one(T)-α)*x.hi )
        C = inv(f_prime(m))

        i += 1
        α /= 2
        i > 20 && error("""Error in guarded_deriv_midpoint:
            the derivative of the function seems too flat""")
    end

    return m, C
end


function K(f::Function, f_prime::Function, x::Interval{T}) where {T}
    m, C = guarded_derivative_midpoint(f, f_prime, x)
    deriv = f_prime(x)
    Kx = m - C*f(m) + (one(T) - C*deriv) * (x - m)
    Kx
end


function krawczyk_refine(f::Function, f_prime::Function, x::Interval{T};
    tolerance=eps(one(T)), debug=false) where {T}

    debug && (print("Entering krawczyk_refine:"); @show x)

    while diam(x) > tolerance  # avoid problem with tiny floating-point numbers if 0 is a root
        Kx = K(f, f_prime, x)
        debug && @show(x, Kx)
        Kx = Kx ∩ x
        Kx == x && break
        x = Kx
    end

    return [Root(x, :unique)]
end




function krawczyk(f::Function, f_prime::Function, x::Interval{T}, level::Int=0;
    tolerance=eps(one(T)), debug=false, maxlevel=30) where {T}

    debug && (print("Entering krawczyk:"); @show(level); @show(x))

    # Maximum level of bisection
    level >= maxlevel && return [Root(x, :unknown)]

    isempty(x) && return Root{typeof(x)}[]  # [(x, :none)]
    Kx = K(f, f_prime, x) ∩ x

    isempty(Kx ∩ x) && return Root{typeof(x)}[]  # [(x, :none)]

    if Kx ⪽ x  # isinterior
        debug && (print("Refining "); @show(x))
        return krawczyk_refine(f, f_prime, Kx, tolerance=tolerance, debug=debug)
    end

    m = mid(x)

    debug && @show(x,m)

    isthin(x) && return [Root(x, :unknown)]

    # bisecting
    roots = vcat(
        krawczyk(f, f_prime, Interval(x.lo, m), level+1,
                 tolerance=tolerance, debug=debug, maxlevel=maxlevel),
        krawczyk(f, f_prime, Interval(m, x.hi), level+1,
                 tolerance=tolerance, debug=debug, maxlevel=maxlevel)
        )


    roots = clean_roots(f, roots)

    return roots
end


# use automatic differentiation if no derivative function given
krawczyk(f::Function,x::Interval{T}; args...) where {T} =
    krawczyk(f, x->D(f,x), x; args...)

# krawczyk for vector of intervals:
krawczyk(f::Function, f_prime::Function, xx::Vector{Interval{T}}; args...) where {T} =
    vcat([krawczyk(f, f_prime, @interval(x); args...) for x in xx]...)

krawczyk(f::Function,  xx::Vector{Interval{T}}, level; args...) where {T} =
    krawczyk(f, x->D(f,x), xx; args...)

krawczyk(f::Function,  xx::Vector{Root{T}}; args...) where {T} =
    krawczyk(f, x->D(f,x), [x.interval for x in xx]; args...)
