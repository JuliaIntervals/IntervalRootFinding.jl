
"""
Make a function ``g: \mathbb{R}^2 \to \mathbb{R}^2`` from
a function ``f: \mathbb{C} \to \mathbb{C}``.

Returns a function taking an `IntervalBox` to an `IntervalBox`.
"""
function realify(f)
    function g(X)
        z = Complex(X[1], X[2])
        zz = f(z)

        return SVector(reim(zz))
    end

    return g
end

"""
    complex_bisection(f, X)

Find complex roots of ``f: \mathbb{C} \to \mathbb{C}``.

Inputs:

- `f`: function that takes ``z \in \mathbb{C}`` and returns another
complex number.

- `X`: An `IntervalBox` specifying the bounds on the real and imaginary parts
of `z`.

"""
function complex_bisection(f, X::IntervalBox, tol=1e-3)

    """
    Make a 2D real version of the complex function `f` suitable for `bisection`,
    i.e. that accepts an `IntervalBox` and returns an `IntervalBox`
    """

    g = realify(f)

    roots = bisection(g, X, tolerance=tol)

    return g, [Complex(root.interval...) for root in roots]
end

"""
    complex_bisection(f, x, y)

Version in which the bounds are specified as two separate `Interval`s.
"""
complex_bisection(f, x::Interval, y::Interval, tol=1e-3) = complex_bisection(f, x × y, tol)
