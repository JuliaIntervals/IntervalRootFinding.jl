# Example of finding roots of a complex function
# f(z) = z^4 + z - 2  (book of Hansen-Walster, chap. 10)

using IntervalRootFinding, IntervalArithmetic

@intervalbox f(x, y) = (z = x + im*y;
                        z2 = z^4 + z - 2;
                        reim(z2)
                        )

roots = bisection(f, (-10..10) × (-10..10) )


