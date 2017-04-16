var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#IntervalRootFinding.jl-1",
    "page": "Home",
    "title": "IntervalRootFinding.jl",
    "category": "section",
    "text": "This Julia package implements guaranteed root-finding methods using interval arithmetic."
},

{
    "location": "index.html#Root-finding-1",
    "page": "Home",
    "title": "Root finding",
    "category": "section",
    "text": "Interval arithmetic not only provides guaranteed numerical calculations; it also makes possible fundamentally new algorithms."
},

{
    "location": "index.html#Interval-Newton-method-1",
    "page": "Home",
    "title": "Interval Newton method",
    "category": "section",
    "text": "One such algorithm is the interval Newton method. This is a version of the standard Newton (or Newton-Raphson) algorithm, an iterative method for finding roots (zeros) of functions. The interval version, however, is fundamentally different from its standard counterpart, in that it can (under the best circumstances) provide rigorous guarantees about the presence or absence and uniqueness of roots of a given function in a given interval, and tells us explicitly when it is unable to provide such a guarantee.The idea of the Newton method is to calculate a root x^ast of a function f [i.e., a value such that f(x^*) = 0] from an initial guess x usingx^* = x - fracf(x)f(xi)for some xi between x and x^*. Since xi is unknown, we can bound it asf(xi) in F(X)where X is a containing interval and F(X) denotes the interval extension of the function f, consisting of applying the same operations as the function f to the interval X.We define an interval Newton operator mathcalN as follows:mathcalN(X) = m(X) - fracF(m(X))F(X)where m(X)  is the midpoint of X converted into an interval.It turns out that mathcalN tells us precisely whether there is a root of f in the interval X: there is no root if mathcalN(X) cap X = emptyset, and there is a unique root if mathcalN(X) subseteq X. There is also an extension to intervals in which the derivative F(X) contains 0, in which case the Newton operator returns a union of two intervals.Iterating the Newton operator on the resulting sets gives a rigorous algorithm that is guaranteed to find all roots of a real function in a given interval (or to inform us if it is unable to do so, for example at a multiple root); see Tucker's book for more details."
},

{
    "location": "index.html#Usage-of-the-interval-Newton-method-1",
    "page": "Home",
    "title": "Usage of the interval Newton method",
    "category": "section",
    "text": "Root-finding routines are in a separate RootFinding submodule of IntervalArithmetic.jl, which must be loaded withjulia> using IntervalArithmetic, IntervalRootFindingThe interval Newton method is implemented for real functions of a single variable as the function newton. For example, we can calculate rigorously the square roots of 2:julia> using ValidatedNumerics\n\njulia> f(x) = x^2 - 2\nf (generic function with 1 method)\n\njulia> newton(f, @interval(-5, 5))\n2-element Array{ValidatedNumerics.Root{Float64},1}:\n Root([-1.4142135623730951, -1.414213562373095], :unique)\n Root([1.414213562373095, 1.4142135623730951], :unique)The function newton  is passed the function and the interval in which to search for roots; it returns an array of Root objects, that contain the interval where a root is found, together with a symbol :unique if there is guaranteed to be a unique root in that interval, or :unknown if the Newton method is unable to make a guarantee, for example, when there is a double root:julia> newton(f, @interval(-5,5))\n6-element Array{ValidatedNumerics.Root{Float64},1}:\n Root([0.9999999968789343, 0.999999997726216], :unknown)\n Root([0.9999999977262161, 0.9999999985734976], :unknown)\n Root([0.9999999987089422, 0.9999999993384274], :unknown)\n Root([0.9999999993384275, 0.9999999999679127], :unknown)\n Root([0.9999999999687099, 1.0000000004524654], :unknown)\n Root([2.0, 2.0], :unique)The Newton method may be applied directly to a vector of known roots, for example to refine them with higher precision:julia> f(x) = x^2 - 2\nf (generic function with 1 method)\n\njulia> roots = newton(f, @interval(-5, 5))\n2-element Array{ValidatedNumerics.Root{Float64},1}:\n Root([-1.4142135623730951, -1.414213562373095], :unique)\n Root([1.414213562373095, 1.4142135623730951], :unique)\n\njulia> setprecision(Interval, 256)\n256\n\njulia> newton(f, roots)\n2-element Array{ValidatedNumerics.Root{Base.MPFR.BigFloat},1}:\n Root([-1.414213562373095048801688724209698078569671875376948073176679737990732478462119, -1.414213562373095048801688724209698078569671875376948073176679737990732478462102]₂₅₆, :unique)\n Root([1.414213562373095048801688724209698078569671875376948073176679737990732478462102, 1.414213562373095048801688724209698078569671875376948073176679737990732478462119]₂₅₆, :unique)\n\njulia> abs(roots2[2].interval.lo - sqrt(big(2)))\n0.000000000000000000000000000000000000000000000000000000000000000000000000000000\n"
},

{
    "location": "index.html#Krawczyk-method-1",
    "page": "Home",
    "title": "Krawczyk method",
    "category": "section",
    "text": "An alternative method is the Krawczyk method, implemented in the function krawczyk, with the same interface as the Newton method:julia> f(x) = x^2 - 2\nf (generic function with 1 method)\n\njulia> krawczyk(f, @interval(-5, 5))\n2-element Array{Root{Float64},1}:\n Root([-1.4142135623730954, -1.4142135623730947], :unique)\n Root([1.4142135623730947, 1.4142135623730954], :unique)\n\njulia> newton(f, @interval(-5, 5))\n2-element Array{Root{Float64},1}:\n Root([-1.4142135623730951, -1.414213562373095], :unique)\n Root([1.414213562373095, 1.4142135623730951], :unique)The Krawczyk method really comes into its own for higher-dimensional functions; this is planned to be implemented in the future."
},

{
    "location": "index.html#find_roots-interface-1",
    "page": "Home",
    "title": "find_roots interface",
    "category": "section",
    "text": "Automatic differentiation is used to calculate the derivative used in the Newton method if the derivative function is not given explicitly as the second argument to newton.An interface find_roots is provided, which does not require an interval to be passed:julia> find_roots(f, -5, 5)\n6-element Array{ValidatedNumerics.Root{Float64},1}:\n Root([0.9999999968789343, 0.999999997726216], :unknown)\n Root([0.9999999977262161, 0.9999999985734976], :unknown)\n Root([0.9999999987089422, 0.9999999993384274], :unknown)\n Root([0.9999999993384275, 0.9999999999679127], :unknown)\n Root([0.9999999999687099, 1.0000000004524654], :unknown)\n Root([1.9999999999999998, 2.0000000000000004], :unique)There is also a version find_roots_midpoint that returns three vectors: the midpoint of each interval; the radius of the interval; and the symbol. This may be useful for someone who just wishes to find roots of a function, without wanting to understand how to manipulate interval objects:julia> find_roots_midpoint(f, -5, 5)\n([-1.4142135623730951,1.414213562373095],[2.220446049250313e-16,4.440892098500626e-16],[:unique,:unique])This uses the function midpoint_radius, that returns the midpoint and radius of a given interval:julia> a = @interval(0.1, 0.2)\n[0.09999999999999999, 0.2]\n\njulia> midpoint_radius(a)\n(0.15,0.05000000000000002)"
},

{
    "location": "api.html#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api.html#IntervalRootFinding.bisect",
    "page": "API",
    "title": "IntervalRootFinding.bisect",
    "category": "Function",
    "text": "bisect(X::IntervalBox, i::Integer, α=0.5)\n\nBisect the IntervalBox in side number i.\n\n\n\n"
},

{
    "location": "api.html#IntervalRootFinding.bisect",
    "page": "API",
    "title": "IntervalRootFinding.bisect",
    "category": "Function",
    "text": "bisect(X::IntervalBox, α=0.5)\n\nBisect the IntervalBox X at position α ∈ [0,1] along its longest side.\n\n\n\n"
},

{
    "location": "api.html#IntervalRootFinding.bisect",
    "page": "API",
    "title": "IntervalRootFinding.bisect",
    "category": "Function",
    "text": "bisect(X::Interval, α=0.5)\n\nSplit the interval X at position α; α=0.5 corresponds to the midpoint. Returns a tuple of the new intervals.\n\n\n\n"
},

{
    "location": "api.html#IntervalRootFinding.newton",
    "page": "API",
    "title": "IntervalRootFinding.newton",
    "category": "Function",
    "text": "newton performs the interval Newton method on the given function f with its optional derivative f_prime and initial interval x. Optional keyword arguments give the tolerance, maxlevel at which to stop subdividing, and a debug boolean argument that prints out diagnostic information.\n\n\n\n"
},

{
    "location": "api.html#IntervalRootFinding.guarded_derivative_midpoint-Tuple{Function,Function,IntervalArithmetic.Interval{T}}",
    "page": "API",
    "title": "IntervalRootFinding.guarded_derivative_midpoint",
    "category": "Method",
    "text": "Returns two intervals, the first being a point within the interval x such that the interval corresponding to the derivative of f there does not contain zero, and the second is the inverse of its derivative\n\n\n\n"
},

{
    "location": "api.html#IntervalRootFinding.guarded_mid-Tuple{Any,IntervalArithmetic.Interval{T}}",
    "page": "API",
    "title": "IntervalRootFinding.guarded_mid",
    "category": "Method",
    "text": "Returns the midpoint of the interval x, slightly shifted in case the midpoint is an exact root\n\n\n\n"
},

{
    "location": "api.html#IntervalRootFinding.newton_refine-Tuple{Function,Function,IntervalArithmetic.Interval{T}}",
    "page": "API",
    "title": "IntervalRootFinding.newton_refine",
    "category": "Method",
    "text": "If a root is known to be inside an interval, newton_refine iterates the interval Newton method until that root is found.\n\n\n\n"
},

{
    "location": "api.html#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": "Modules = [IntervalRootFinding]\nOrder   = [:type, :function]"
},

]}
