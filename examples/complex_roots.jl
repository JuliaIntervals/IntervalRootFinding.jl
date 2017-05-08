# Example of finding roots of a complex function
# f(z) = z^4 + z - 2  (book of Hansen-Walster, chap. 10)

using IntervalRootFinding, IntervalArithmetic

doc"""
    complex_bisection(f, X)

Find complex roots of $f: \mathbb{C} \to \mathbb{C}$.

Inputs:

- `f`: function $f: \mathbb{C} \to \mathbb{C}$, i.e. a function that accepts a complex number and returns another complex number.

- `X`: An `IntervalBox` specifying the bounds on the real and imaginary parts
of `z`.

"""
function complex_bisection(f, X::IntervalBox)

    """
    Make a 2D real version of the complex function `f` suitable for `bisection`,
    i.e. that accepts an `IntervalBox` and returns an `IntervalBox`
    """
    function g(X::IntervalBox)
        z = Complex(X...)
        zz = f(z)

        return IntervalBox(reim(zz))
    end

    roots = bisection(g, X)

    return g, [Complex(root.interval...) for root in roots]
end

"""
    complex_bisection(f, x, y)

Version in which the bounds are specified as two separate `Interval`s.
"""
complex_bisection(f, x::Interval, y::Interval) = complex_bisection(f, x × y)

f(z) = z^4 + z - 2

L = 10
g, roots = complex_bisection(f, -L..L, -L..L)

println(roots)
