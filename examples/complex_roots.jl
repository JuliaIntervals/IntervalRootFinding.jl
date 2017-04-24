# Example of finding roots of a complex function
# f(z) = z^4 + z - 2  (book of Hansen-Walster, chap. 10)

using IntervalRootFinding, IntervalArithmetic

# Find complex roots by bisection
function complex_bisection(f, X::IntervalBox)

    function g(X::IntervalBox)
        x, y = X
        z = x + im*y;
        zz = f(z)
        x, y = reim(zz)
        return IntervalBox(x, y)
    end

    roots = bisection(g, X)

    return [root.interval[1] + im*root.interval[2] for root in roots]
end

f(z) = z^4 + z - 2

complex_bisection(f, (-100..100) × (-100..100) )

