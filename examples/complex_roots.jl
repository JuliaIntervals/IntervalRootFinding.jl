# Example of finding roots of a complex function
# f(z) = z^4 + z - 2  (book of Hansen-Walster, chap. 10)

using IntervalRootFinding, IntervalArithmetic

"""
    complex_bisection(f, X)

Find complex roots of $f: \mathbb{C} \to \mathbb{C}$.

Inputs:

- `f`: function that takes $z \in \mathbb{C}$ and returns another
complex number.

- `X`: An `IntervalBox` specifying the bounds on the real and imaginary parts
of `z`.

"""
function complex_bisection(f, X::IntervalBox)

    """
    Make a 2D real version of the complex function `f` suitable for `bisection`,
    i.e. that accepts an `IntervalBox` and returns an `IntervalBox`
    """
    function g(X::IntervalBox)
        x, y = X
        z = x + im*y;
        zz = f(z)
        x, y = reim(zz)
        return IntervalBox(x, y)
    end

    roots = bisection(g, X)

    return g, [root.interval[1] + im*root.interval[2] for root in roots]
end

"""
    complex_bisection(f, x, y)

Version in which the bounds are specified as two separate `Interval`s.
"""
complex_bisection(f, x::Interval, y::Interval) = complex_bisection(f, x × y)

f(z) = z^4 + z - 2

L = 10
g, roots = complex_bisection(f, -L..L, -L..L)

