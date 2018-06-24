# This file is part of the ValidatedNumerics.jl package; MIT licensed

# Newton method



doc"""If a root is known to be inside an interval,
`newton_refine` iterates the interval Newton method until that root is found."""
function newton_refine{N,T}(f::Function, f_prime::Function, X::Union{Interval{T}, IntervalBox{N,T}};
                          tolerance=eps(T), debug=false)

    debug && (print("Entering newton_refine:"); @show X)

    while diam(X) > tolerance  # avoid problem with tiny floating-point numbers if 0 is a root
        deriv = f_prime(X)
        NX = 𝒩(f, X, deriv)

        debug && @show(X, NX)

        NX = NX ∩ X
        NX == X && break

        X = NX
    end

    debug && @show "Refined root", X

    return [Root(X, :unique)]
end



doc"""If a root is known to be inside an interval,
`newton_refine` iterates the interval Newton method until that root is found."""
function newton_refine{T}(f::Function, f_prime::Function, X::Interval{T};
                          tolerance=eps(T), debug=false)

    debug && (print("Entering newton_refine:"); @show X)

    while diam(X) > tolerance  # avoid problem with tiny floating-point numbers if 0 is a root
        deriv = f_prime(X)
        NX = 𝒩(f, X, deriv)

        debug && @show(X, NX)

        NX = NX ∩ X
        NX == X && break

        X = NX
    end

    debug && @show "Refined root", X

    return [Root(X, :unique)]
end







doc"""`newton` performs the interval Newton method on the given function `f`
with its optional derivative `f_prime` and initial interval `x`.
Optional keyword arguments give the `tolerance`, `maxlevel` at which to stop
subdividing, and a `debug` boolean argument that prints out diagnostic information."""

function newton{T}(f::Function, f_prime::Function, x::Interval{T}, level::Int=0;
                   tolerance=eps(T), debug=false, maxlevel=30)

    debug && (print("Entering newton:"); @show(level); @show(x))

    isempty(x) && return Root{typeof(x)}[]  #[(x, :none)]

    z = zero(x.lo)
    tolerance = abs(tolerance)

    # Tolerance reached or maximum level of bisection
    (diam(x) < tolerance || level >= maxlevel) && return [Root(x, :unknown)]

    deriv = f_prime(x)
    debug && @show(deriv)

    if !(z in deriv)
        Nx = 𝒩(f, x, deriv)
        debug && @show(Nx, Nx ⊆ x, Nx ∩ x)

        isempty(Nx ∩ x) && return Root{typeof(x)}[]

        if Nx ⊆ x
            debug && (print("Refining "); @show(x))
            return newton_refine(f, f_prime, Nx, tolerance=tolerance, debug=debug)
        end

        m = mid(x)

        debug && @show(x,m)

        # bisect:
        roots = vcat(
            newton(f, f_prime, Interval(x.lo, m), level+1,
                   tolerance=tolerance, debug=debug, maxlevel=maxlevel),
            newton(f, f_prime, Interval(m, x.hi), level+1,
                   tolerance=tolerance, debug=debug, maxlevel=maxlevel)
            # note the nextfloat here to prevent repetition
            )

        debug && @show(roots)

    else  # 0 in deriv; this does extended interval division by hand

        y1 = Interval(deriv.lo, -z)
        y2 = Interval(z, deriv.hi)

        if debug
            @show (y1, y2)
            @show N(f, x, y1)
            @show N(f, x, y2)
        end

        y1 = 𝒩(f, x, y1) ∩ x
        y2 = 𝒩(f, x, y2) ∩ x

        debug && @show(y1, y2)

        roots = vcat(
            newton(f, f_prime, y1, level+1,
                   tolerance=tolerance, debug=debug, maxlevel=maxlevel),
            newton(f, f_prime, y2, level+1,
                   tolerance=tolerance, debug=debug, maxlevel=maxlevel)
            )

        debug && @show(roots)

    end


    roots = clean_roots(f, roots)

    return roots

end


# use automatic differentiation if no derivative function given:
newton{T}(f::Function, x::Interval{T};  args...) =
    newton(f, x->D(f,x), x; args...)

# newton for vector of intervals:
newton{T}(f::Function, f_prime::Function, xx::Vector{Interval{T}}; args...) =
    vcat([newton(f, f_prime, @interval(x); args...) for x in xx]...)

newton{T}(f::Function,  xx::Vector{Interval{T}}, level; args...) =
    newton(f, x->D(f,x), xx, 0, args...)

newton{T}(f::Function,  xx::Vector{Root{T}}; args...) =
    newton(f, x->D(f,x), [x.interval for x in xx], args...)
