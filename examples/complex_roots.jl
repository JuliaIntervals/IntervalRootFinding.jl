# Example of finding roots of a complex function
# f(z) = z^4 + z - 2  (book of Hansen-Walster, chap. 10)

using IntervalRootFinding, IntervalArithmetic



f(z) = z^4 + z - 2

L = 10
g, roots = complex_bisection(f, -L..L, -L..L)
