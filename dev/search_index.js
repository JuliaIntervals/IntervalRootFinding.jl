var documenterSearchIndex = {"docs":
[{"location":"roots/#roots-interface","page":"roots interface","title":"roots interface","text":"","category":"section"},{"location":"roots/#Methods","page":"roots interface","title":"Methods","text":"","category":"section"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"Three root-finding contractors currently available through the roots interface are the following:","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"Newton (default);\nKrawczyk;\nBisection","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"Both the Newton and Krawczyk methods can determine if a root is unique in an interval, at the cost of requiring that the function is differentiable. The bisection method has no such requirement, but can never guarantee the existence or uniqueness of a root.","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"The method used is given using the contractor keyword argument:","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"julia> roots(log, -2..2 ; contractor = Newton)\n1-element Vector{Root{Interval{Float64}}}:\n Root([0.999999, 1.00001]_com, :unique)\n\njulia> roots(log, -2..2 ; contractor = Krawczyk)\n1-element Vector{Root{Interval{Float64}}}:\n Root([0.999999, 1.00001]_com, :unique)\n\njulia> roots(log, -2..2 ; contractor = Bisection)\n1-element Vector{Root{Interval{Float64}}}:\n Root([0.999999, 1.00001]_com, :unknown)","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"Note that as shown in the example, the log function does not complain about being given an interval going outside of its domain. While this may be surprising, this is the expected behavior and no root will ever be found outside the domain of a function.","category":"page"},{"location":"roots/#Explicit-derivatives","page":"roots interface","title":"Explicit derivatives","text":"","category":"section"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"Newton and Krawczyk methods require the function to be differentiable, but the derivative is usually computed automatically using forward-mode automatic differentiation, provided by the ForwardDiff.jl package. It is however possible to provide the derivative explicitly using the derivative keyword argument:","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"julia> roots(log, -2..2 ; contractor = Newton, derivative = x -> 1/x)\n1-element Vector{Root{Interval{Float64}}}:\n Root([0.999999, 1.00001]_com_NG, :unique)\n\njulia> roots(log, -2..2 ; contractor = Krawczyk, derivative = x -> 1/x)\n1-element Vector{Root{Interval{Float64}}}:\n Root([0.999999, 1.00001]_com_NG, :unique)","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"When providing the derivative explicitly, the computation is expected to be slightly faster, but the precision of the result is unlikely to be affected.","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"julia> using BenchmarkTools\n\njulia> @btime roots(log, -2..2 ; derivative = x -> 1/x)\n  4.814 μs (53 allocations: 3.16 KiB)\n1-element Vector{Root{Interval{Float64}}}:\n Root([0.999999, 1.00001]_com_NG, :unique)\n\njulia> @btime roots(log, -2..2)\n  5.767 μs (53 allocations: 3.16 KiB)\n1-element Vector{Root{Interval{Float64}}}:\n Root([0.999999, 1.00001]_com, :unique)","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"This may be useful in some special cases where ForwardDiff.jl is unable to compute the derivative of a function. Examples are complex functions and functions whose interval extension must be manually defined (e.g. special functions like zeta).","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"In dimension greater than one, the derivative be given as a function returning the Jacobi matrix:","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"julia> function f( (x, y) )\n           return [sin(x), cos(y)]\n       end\nf (generic function with 1 method)\n\njulia> function df( (x, y) )\n           return [cos(x) 0 ; 0 -sin(y)]\n       end\n\njulia> roots(f, [-3..3, -3..3] ; derivative = df)\n2-element Vector{Root{Vector{Interval{Float64}}}}:\n Root(Interval{Float64}[[-1.24409e-21, 1.0588e-22]_com_NG, [-1.5708, -1.57079]_com_NG], :unique)\n Root(Interval{Float64}[[-1.24409e-21, 1.0588e-22]_com_NG, [1.57079, 1.5708]_com_NG], :unique)","category":"page"},{"location":"roots/#Tolerance","page":"roots interface","title":"Tolerance","text":"","category":"section"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"An absolute tolerance for the search may be specified as the abstol keyword argument. Currently a method must first be provided in order to be able to choose the tolerance.","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"julia> g(x) = sin(exp(x))\ng (generic function with 1 method)\n        \njulia> roots(g, 0..2)\n2-element Vector{Root{Interval{Float64}}}:\n Root([1.14472, 1.14474]_com, :unique)\n Root([1.83787, 1.83788]_com, :unique)\n\njulia> roots(g, 0..2 ; abstol = 1e-1)\n2-element Vector{Root{Interval{Float64}}}:\n Root([1.14173, 1.15244]_com, :unique)\n Root([1.78757, 1.84273]_com, :unique)","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"A lower tolerance may greatly reduce the computation time, at the cost of an increased number of returned roots having :unknown status:","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"julia> h(x) = cos(x) * sin(1 / x)\nh (generic function with 1 method)\n\njulia> @btime roots(h, 0.05..1)\n  79.500 μs (301 allocations: 13.97 KiB)\n6-element Vector{Root{Interval{Float64}}}:\n Root([0.0530516, 0.0530517]_com_NG, :unique)\n Root([0.0636619, 0.063662]_com_NG, :unique)\n Root([0.0795774, 0.0795775]_com_NG, :unique)\n Root([0.106103, 0.106104]_com_NG, :unique)\n Root([0.159154, 0.159156]_com_NG, :unique)\n Root([0.318309, 0.31831]_com_NG, :unique)\n\njulia> @btime roots(h, 0.05..1 ; abstol = 1e-2)\n  48.500 μs (265 allocations: 11.61 KiB)\n6-element Vector{Root{Interval{Float64}}}:\n Root([0.0514445, 0.0531087]_com_NG, :unique)\n Root([0.0570253, 0.0641615]_com, :unknown)\n Root([0.0785458, 0.0797797]_com_NG, :unknown)\n Root([0.104754, 0.107542]_com_NG, :unknown)\n Root([0.157236, 0.165989]_com_NG, :unknown)\n Root([0.31716, 0.319318]_com_NG, :unique)\n\njulia> @btime roots(h, 0.05..1 ; abstol = 1e-1)\n  17.300 μs (114 allocations: 5.16 KiB)\n3-element Vector{Root{Interval{Float64}}}:\n Root([0.04999999, 0.107542]_com, :unknown)\n Root([0.107541, 0.165989]_com, :unknown)\n Root([0.283803, 0.356099]_com_NG, :unknown)","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"The last example shows a case where the tolerance was too large to be able to isolate the roots in distinct regions.","category":"page"},{"location":"roots/","page":"roots interface","title":"roots interface","text":"warning: Warning\nFor a root x of some function, if the absolute tolerance is smaller than eps(x) i.e. if tol + x == x, roots may never be able to converge to the required tolerance and the function may get stuck in an infinite loop.","category":"page"},{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/#Main-interface","page":"API","title":"Main interface","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [IntervalRootFinding]\nPages   = [\"roots.jl\"]","category":"page"},{"location":"api/#IntervalRootFinding.RootProblem-Tuple{Any, Any}","page":"API","title":"IntervalRootFinding.RootProblem","text":"RootProblem(f::Function, search_region ; kwargs...)\n\nSetup a RootProblem for searching the roots (also known as zeros) of the function f in the given search region.\n\nsearch_region must be either an Interval if f is scalar or a vector of Intervals if f is vector-valued.\n\nThe returned RootProblem is an iterator that give access to the internal state of the search during the iteration, allowing to add callbacks and logging to the search.\n\nParameters\n\ncontractor: Contractor used to determine the status of a region.   Must be either Newton, Krawczyk, or Bisection. Bisection do not require   to compute the derivative of the function, but can also never guarantee the   existence of a root. Default: Newton.\nderivative: Explicit derivative of the function (or its jacobian for   vector-valued functions) used by the Newton and Krawczyk contractors.   Default: nothing (the derivative is computed automatically using ForwardDiff.jl).\nsearch_order: Order in which the sub-regions are searched.   BreadthFirst (visit the largest regions first) and DepthFirst   (visit the smallest regions first) are supported. Default: BreadthFirst.\nabstol: Absolute tolerance. The search is stopped when all dimension   of the remaining regions are smaller than abstol. Default: 1e-7.\nwhere_bisect: Value used to bisect the region. It is used to avoid   bisecting exactly on zero when starting with symmetrical regions,   often leading to having a solution directly on the boundary of a region,   which prevent the contractor to prove it's unicity. Default: 127/256.\n\nreltol and max_iteration are currently ignored.\n\n\n\n\n\n","category":"method"},{"location":"api/#IntervalRootFinding.roots-Tuple{Any, Any}","page":"API","title":"IntervalRootFinding.roots","text":"roots(f::Function, search_region ; kwargs...)\n\nReturn the roots (also known as zeros) of the function f contained in the given search region (either an Interval for a scalar function or vector of Intervals for a vector valued function), together with a status.\n\nThe status of the returned regions can be either of\n\n:unique: the region provably contains exactly one root of f\n:unknown: the region may contain any number of roots (potentially none)\n\nThe parts of the search region that are not covered by any of the returned roots are guaranteed to contain no root of the function.\n\nFor information about the optional search parameters, see RootProblem.\n\n\n\n\n\n","category":"method"},{"location":"api/#Root-object","page":"API","title":"Root object","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [IntervalRootFinding]\nPages   = [\"root_object.jl\"]","category":"page"},{"location":"api/#IntervalRootFinding.Root","page":"API","title":"IntervalRootFinding.Root","text":"Root\n\nObject representing a possible root inside a given region. The field status is either :unknown or :unique. If status is :unique then we know that there is a unique root of the function in question inside the given region.\n\nInternally the status may also be :empty for region guaranteed to contain no root, however such Roots are discarded by default and thus never returned by the roots function.\n\nFields\n\nregion: a region (either Interval or SVector of interval     representing an interval box) searched for roots.\nstatus: the status of the region, valid values are :empty, unknown     and :unique.\n\n\n\n\n\n","category":"type"},{"location":"api/#IntervalRootFinding.isunique-Union{Tuple{Root{T}}, Tuple{T}} where T","page":"API","title":"IntervalRootFinding.isunique","text":"isunique(rt)\n\nReturn whether a Root is unique.\n\n\n\n\n\n","category":"method"},{"location":"api/#IntervalRootFinding.root_region-Tuple{Root}","page":"API","title":"IntervalRootFinding.root_region","text":"root_region(rt)\n\nReturn the region associated to a Root.\n\n\n\n\n\n","category":"method"},{"location":"api/#IntervalRootFinding.root_status-Tuple{Root}","page":"API","title":"IntervalRootFinding.root_status","text":"root_status(rt)\n\nReturn the status of a Root.\n\n\n\n\n\n","category":"method"},{"location":"api/#Contractors","page":"API","title":"Contractors","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [IntervalRootFinding]\nPages   = [\"contractors.jl\"]","category":"page"},{"location":"api/#IntervalRootFinding.AbstractContractor","page":"API","title":"IntervalRootFinding.AbstractContractor","text":"AbstractContractor\n\nAbstract type for contractors.\n\n\n\n\n\n","category":"type"},{"location":"api/#IntervalRootFinding.refine-Union{Tuple{C}, Tuple{RootProblem{C}, Root}} where C","page":"API","title":"IntervalRootFinding.refine","text":"refine(root_problem::RootProblem, X::Root)\n\nRefine a root.\n\n\n\n\n\n","category":"method"},{"location":"api/#Branch-and-bound-search-interface","page":"API","title":"Branch-and-bound search interface","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [IntervalRootFinding]\nPages   = [\"branch_and_bound.jl\"]","category":"page"},{"location":"api/#Others","page":"API","title":"Others","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [IntervalRootFinding]\nPages   = [\"complex.jl\", \"linear_eq.jl\", \"newton1d.jl\", \"quadratic.jl\", \"slopes.jl\"]","category":"page"},{"location":"api/#IntervalRootFinding.realify_derivative-Tuple{Any}","page":"API","title":"IntervalRootFinding.realify_derivative","text":"Takes the derivative of a complex function and returns the real jacobian that implements it.\n\n\n\n\n\n","category":"method"},{"location":"api/#IntervalRootFinding.gauss_elimination_interval!-Tuple{AbstractArray, AbstractMatrix, AbstractArray}","page":"API","title":"IntervalRootFinding.gauss_elimination_interval!","text":"Solves the system of linear equations using Gaussian Elimination. Preconditioning is used when the precondition keyword argument is true.\n\nREF: Luc Jaulin et al., Applied Interval Analysis, pg. 72\n\n\n\n\n\n","category":"method"},{"location":"api/#IntervalRootFinding.gauss_elimination_interval1!-Tuple{AbstractArray, AbstractMatrix, AbstractArray}","page":"API","title":"IntervalRootFinding.gauss_elimination_interval1!","text":"Using Base.`\n\n\n\n\n\n","category":"method"},{"location":"api/#IntervalRootFinding.gauss_seidel_interval!-Tuple{AbstractArray, AbstractMatrix, AbstractArray}","page":"API","title":"IntervalRootFinding.gauss_seidel_interval!","text":"Iteratively solves the system of interval linear equations and returns the solution set. Uses the Gauss-Seidel method (Hansen-Sengupta version) to solve the system. Keyword precondition to turn preconditioning off. Eldon Hansen and G. William Walster : Global Optimization Using Interval Analysis - Chapter 5 - Page 115\n\n\n\n\n\n","category":"method"},{"location":"api/#IntervalRootFinding.preconditioner-Tuple{AbstractMatrix, AbstractArray}","page":"API","title":"IntervalRootFinding.preconditioner","text":"Preconditions the matrix A and b with the inverse of mid(A)\n\n\n\n\n\n","category":"method"},{"location":"api/#IntervalRootFinding.newton1d-Union{Tuple{T}, Tuple{Function, Function, Interval{T}}} where T","page":"API","title":"IntervalRootFinding.newton1d","text":"newton1d performs the interval Newton method on the given function f with its derivative f′ and initial interval x. Optional keyword arguments give the tolerances reltol and abstol. reltol is the tolerance on the relative error whereas abstol is the tolerance on |f(X)|, and a debug boolean argument that prints out diagnostic information.\n\n\n\n\n\n","category":"method"},{"location":"api/#IntervalRootFinding.newton1d-Union{Tuple{T}, Tuple{Function, Interval{T}}} where T","page":"API","title":"IntervalRootFinding.newton1d","text":"newton1d performs the interval Newton method on the given function f and initial interval x. Optional keyword arguments give the tolerances reltol and abstol. reltol is the tolerance on the relative error whereas abstol is the tolerance on |f(X)|, and a debug boolean argument that prints out diagnostic information.\n\n\n\n\n\n","category":"method"},{"location":"api/#IntervalRootFinding.quadratic_helper!-Union{Tuple{T}, Tuple{Interval{T}, Interval{T}, Interval{T}, Array{Interval{T}}}} where T","page":"API","title":"IntervalRootFinding.quadratic_helper!","text":"Helper function for quadratic_interval that computes roots of a real quadratic using interval arithmetic to bound rounding errors.\n\n\n\n\n\n","category":"method"},{"location":"api/#IntervalRootFinding.quadratic_roots-Union{Tuple{T}, Tuple{Interval{T}, Interval{T}, Interval{T}}} where T","page":"API","title":"IntervalRootFinding.quadratic_roots","text":"Function to solve a quadratic equation where the coefficients are intervals. Returns an array of intervals of the roots. Arguments a, b and c are interval coefficients of x², x and 1 respectively. The interval case differs from the non-interval case in that there might be three disjoint interval roots. In the third case, one interval root extends to −∞ and another extends to +∞. This algorithm finds the set of points where F.lo(x) ≥ 0 and the set of points where F.hi(x) ≤ 0 and takes the intersection of these two sets. Eldon Hansen and G. William Walster : Global Optimization Using Interval Analysis - Chapter 8\n\n\n\n\n\n","category":"method"},{"location":"decorations/#Decorations-and-guarantees","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"","category":"section"},{"location":"decorations/","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"When you define a simple function, you will notice that the roots usually come out with the flag com and NG.","category":"page"},{"location":"decorations/","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"julia> roots(x -> x^2 - 0.1, interval(-5, 5))\n2-element Vector{Root{Interval{Float64}}}:\n Root([-0.316228, -0.316227]_com_NG, :unique)\n Root([0.316227, 0.316228]_com_NG, :unique)","category":"page"},{"location":"decorations/","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"In this case, com is the decoration of the interval, and NG its guarantee. See the documentation of IntervalArithmetic.jl for a detailed description, here we will go on what they mean in the context of finding roots of a function.","category":"page"},{"location":"decorations/#Decorations","page":"Decorations and guarantees","title":"Decorations","text":"","category":"section"},{"location":"decorations/","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"The decoration tells us whether all operations that affected the interval were well-behaved. The two crucial ones for us are def (the operations were not continuous) and trv (something horribly wrong happened).","category":"page"},{"location":"decorations/","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"In both cases, it means that the root can not be trusted, as the hypothesis of the theory we use are not fulfilled.","category":"page"},{"location":"decorations/#Guarantee","page":"Decorations and guarantees","title":"Guarantee","text":"","category":"section"},{"location":"decorations/","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"The NG flags means that non-interval have been mixed with intervals. Since we can not track the source of the non-interval numbers, we can not guarantee that they are correct. For example, 0.1 is famously not parsed as 0.1, as 0.1 can not be represented exactly as a binary number (just like 1/3 can not be represented exactly as a decimal number).","category":"page"},{"location":"decorations/","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"julia> big(0.1)\n0.1000000000000000055511151231257827021181583404541015625","category":"page"},{"location":"decorations/","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"If you want the number that you are inputting to be trusted as is, you can use the @exact macro from IntervalArithmetic.jl, and the NG flag will be avoided in most cases.","category":"page"},{"location":"decorations/","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"julia> @exact f(x) = x^2 - 0.1\nf (generic function with 1 method)\n\njulia> roots(f, interval(-5, 5))\n2-element Vector{Root{Interval{Float64}}}:\n Root([-0.316228, -0.316227]_com, :unique)\n Root([0.316227, 0.316228]_com, :unique)","category":"page"},{"location":"decorations/","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"The NG flag can still appear if other computations, from a library using non-interval \"magic\" numbers for example, thus indicating that some non-trusted numbers have been used in the computation.","category":"page"},{"location":"decorations/","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"Moreover, the macro work in such a way that you can still use the defined function with numbers and get floating point results.","category":"page"},{"location":"decorations/","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"julia> f(0.2)\n-0.06","category":"page"},{"location":"decorations/#Trust","page":"Decorations and guarantees","title":"Trust","text":"","category":"section"},{"location":"decorations/","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"We are doing our best to really give validated and guaranteed results. However, in the case that you may make your house explode based on a result returned by this package, we would like to remind you that you should not trust the package beyond the promises of the license:","category":"page"},{"location":"decorations/","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED","category":"page"},{"location":"decorations/","page":"Decorations and guarantees","title":"Decorations and guarantees","text":"All bug reports and pull requests to fix them are however more than welcome.","category":"page"},{"location":"#IntervalRootFinding.jl","page":"Home","title":"IntervalRootFinding.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package provides guaranteed methods for finding roots of functions f mathbbR^n to mathbbR^n with n ge 1, i.e. vectors (or scalars, for n=1) mathbbx for which f(mathbbx) = mathbb0. In principle, it guarantees to find all roots inside a given box in mathbbR^n, or report subboxes for which it is unable to provide guarantees.","category":"page"},{"location":"","page":"Home","title":"Home","text":"To do so, it uses methods from interval analysis, using interval arithmetic from the IntervalArithmetic.jl package by the same authors.","category":"page"},{"location":"","page":"Home","title":"Home","text":"warning: Warning\nWhile this package aimed at providing guaranteed results and despite our best efforts and test suite, some bugs may remain and there are several known issues with corner cases. Please look at the issue tracker and report there any odd and/or incorrect behavior.","category":"page"},{"location":"#Basic-1D-example","page":"Home","title":"Basic 1D example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To begin, we need a standard Julia function and an interval in which to search roots of that function. Intervals use the Interval type provided by the IntervalArithmetic.jl package Intervals are generally constructed using the .. syntax (from the IntervalArithmetic.Symbols submodule), a..b representing the closed interval a b.","category":"page"},{"location":"","page":"Home","title":"Home","text":"When provided with this information, the roots function will return a vector of all roots of the function in the given interval.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Example:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using IntervalArithmetic, IntervalArithmetic.Symbols, IntervalRootFinding\n\njulia> rts = roots(x -> x^2 - 2x, 0..10)\n2-element Vector{Root{Interval{Float64}}}:\n Root([0.0, 3.73849e-08]_com_NG, :unknown)\n Root([1.999999, 2.00001]_com_NG, :unique)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The roots are returned as Root objects, containing an interval and the status of that interval, represented as a Symbol. There are two possible types of root status, as shown in the example:","category":"page"},{"location":"","page":"Home","title":"Home","text":":unique: the given interval contains exactly one root of the function,\n:unknown: the given interval may or may not contain one or more roots; the algorithm used was unable to come to a conclusion.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The second status is still informative, since all regions of the original search interval not contained in any of the returned root intervals is guaranteed not to contain any root of the function. In the above example, we know that the function has no root in the interval 21 10, for example.","category":"page"},{"location":"","page":"Home","title":"Home","text":"There are several known situations where the uniqueness (and existence) of a solution cannot be determined by the interval algorithms used in the package:","category":"page"},{"location":"","page":"Home","title":"Home","text":"If the solution is on the boundary of the interval (as in the previous example);\nIf the derivative of the solution is zero at the solution.","category":"page"},{"location":"","page":"Home","title":"Home","text":"In particular, the second condition means that multiple roots cannot be proven to be unique. For example:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> g(x) = (x^2 - 2)^2 * (x^2 - 3)\ng (generic function with 1 method)\n\njulia> roots(g, -10..10)\n4-element Vector{Root{Interval{Float64}}}:\n Root([-1.73206, -1.73205]_com_NG, :unique)\n Root([-1.41422, -1.41421]_com, :unknown)  \n Root([1.41421, 1.41422]_com, :unknown)    \n Root([1.73205, 1.73206]_com_NG, :unique)  ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here we see that the two double roots are reported as being possible roots without guarantee and the simple roots have been proved to be unique.","category":"page"},{"location":"#Basic-multi-dimensional-example","page":"Home","title":"Basic multi-dimensional example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For dimensions n  1, the function passed to roots must take an array as argument and return an array. The initial search region is an array of interval.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here we give a 3D example:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> function g( (x1, x2, x3) )\n          return [\n              x1^2 + x2^2 + x3^2 - 1,\n              x1^2 + x3^2 - 0.25,\n              x1^2 + x2^2 - 4x3\n          ]\n       end\ng (generic function with 1 method)\n\njulia> X = -5..5\n[-5.0, 5.0]_com\n\njulia> rts = roots(g, [X, X, X])\n4-element Vector{Root{Vector{Interval{Float64}}}}:\n Root(Interval{Float64}[[-0.440763, -0.440762]_com_NG, [-0.866026, -0.866025]_com_NG, [0.236067, 0.236069]_com_NG], :unique)\n Root(Interval{Float64}[[-0.440763, -0.440762]_com_NG, [0.866025, 0.866026]_com_NG, [0.236067, 0.236069]_com_NG], :unique)\n Root(Interval{Float64}[[0.440762, 0.440763]_com_NG, [-0.866026, -0.866025]_com_NG, [0.236067, 0.236069]_com_NG], :unique)\n Root(Interval{Float64}[[0.440762, 0.440763]_com_NG, [0.866025, 0.866026]_com_NG, [0.236067, 0.236069]_com_NG], :unique)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Thus, the system admits four unique roots in the box -5 5^3.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Moreover, the package is compatible with StaticArrays.jl. Usage of static arrays is recommended to increase performance.","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using StaticArrays\n\njulia> h((x, y)) = SVector(x^2 - 4, y^2 - 16)\nh (generic function with 1 method)\n\njulia> roots(h, SVector(X, X))\n4-element Vector{Root{SVector{2, Interval{Float64}}}}:\n Root(Interval{Float64}[[-2.00001, -1.999999]_com_NG, [-4.00001, -3.999999]_com_NG], :unique)\n Root(Interval{Float64}[[-2.00001, -1.999999]_com_NG, [3.999999, 4.00001]_com_NG], :unique)\n Root(Interval{Float64}[[1.999999, 2.00001]_com_NG, [-4.00001, -3.999999]_com_NG], :unique)\n Root(Interval{Float64}[[1.999999, 2.00001]_com_NG, [3.999999, 4.00001]_com_NG], :unique)","category":"page"},{"location":"#Stationary-points","page":"Home","title":"Stationary points","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Stationary points of a function fmathbbR^n to mathbbR may be found as zeros of the gradient of f. The gradient can be computed using ForwardDiff.jl:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using ForwardDiff: gradient\n\njulia> f( (x, y) ) = sin(x) * sin(y)\nf (generic function with 1 method)\n\njulia> ∇f(x) = gradient(f, x)  # gradient operator from the package\n(::#53) (generic function with 1 method)\n\njulia> rts = roots(∇f, IntervalBox(-5..6, 2), Newton, 1e-5)\n25-element Vector{Root{Vector{Interval{Float64}}}}:\n Root(Interval{Float64}[[-4.7124, -4.71238]_com, [-4.7124, -4.71238]_com], :unique)\n Root(Interval{Float64}[[-3.1416, -3.14159]_com, [-3.1416, -3.14159]_com], :unique)\n ⋮\n [output skipped for brevity]","category":"page"},{"location":"","page":"Home","title":"Home","text":"Now let's find the midpoints and plot them:","category":"page"},{"location":"","page":"Home","title":"Home","text":"midpoints = [mid.(root_region(rt)) for rt in rts]\n\nxs = first.(midpoints)\nys = last.(midpoints)\n\nusing Plots; plotlyjs()\n\nsurface(-5:0.1:6, -6:0.1:6, (x,y) -> f([x, y]))\nscatter!(xs, ys, f.(midpoints))","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: stationary points)","category":"page"},{"location":"biblio/#Bibliography","page":"Bibliography","title":"Bibliography","text":"","category":"section"},{"location":"biblio/","page":"Bibliography","title":"Bibliography","text":"Applied Interval Analysis, L. Jaulin, M. Kieffer, O. Didrit & E. Walter. Springer (2001)\nGlobal Optimization Using Interval Analysis: Revised And Expanded, E. Hansen & G.W. Walster. CRC Press (2003)\nIntroduction to Interval Analysis, R.E. Moore, R.B. Kearfott & M.J. Cloud. SIAM (2009)\nValidated Numerics: A Short Introduction to Rigorous Computations, W. Tucker. Princeton University Press (2010)","category":"page"},{"location":"internals/#Internals","page":"Internals","title":"Internals","text":"","category":"section"},{"location":"internals/","page":"Internals","title":"Internals","text":"This section describes some of the internal mechanism of the package and several ways to use them to customize a search.","category":"page"},{"location":"internals/#Branch-and-bound","page":"Internals","title":"Branch and bound","text":"","category":"section"},{"location":"internals/","page":"Internals","title":"Internals","text":"When roots is called, it performs a branch-and-bound search, that is, it iteratively looks at a region X and for each region tries to determine if it contains a root. It then does the following:","category":"page"},{"location":"internals/","page":"Internals","title":"Internals","text":"If X is proven to contain no root, it discards it.\nIf X is proven to contain exactly one root, it tries to get the best possible bounds for the root and store the resulting region in the list of Roots to output with a :unique status.\nIf the test is inconclusive and the size of X is smaller than the tolerance, it stores it in the list of Roots to output with :unknown status.\nIf the test is inconclusive and the size of X is larger than the tolerance, it bisects X and then processes each resulting half.","category":"page"},{"location":"internals/","page":"Internals","title":"Internals","text":"At some point all regions will either have a determined status or be smaller than the tolerance, and the algorithm will halt and return all stored roots.","category":"page"},{"location":"internals/#Tree-representation","page":"Internals","title":"Tree representation","text":"","category":"section"},{"location":"internals/","page":"Internals","title":"Internals","text":"A branch-and-bound search can be naturally represented as a binary tree: each leaf contains a region and its status and each node represents a bisection. If the tree is known, the topology of the whole search can be reconstructed. There is however a point not determined by the tree.","category":"page"},{"location":"internals/","page":"Internals","title":"Internals","text":"The branch and bound algorithm used by the roots function builds this tree and at the end collect all leaves containing a region with status either :unknown or :unique. We see below how to access the tree.","category":"page"},{"location":"internals/#Search-strategy","page":"Internals","title":"Search strategy","text":"","category":"section"},{"location":"internals/","page":"Internals","title":"Internals","text":"While the tree representation is sufficient to know the topology of the search, it does not determine the order in which leaves are processed during the search. This has an influence on how the tree is created and the amount of memory used, but will not change the resulting tree, unless some limitations on the number of iterations or leaves are enforced.","category":"page"},{"location":"internals/","page":"Internals","title":"Internals","text":"note: Note\nNo such limitations are currently implemented, but they are planned. They will allow to deal, for example, with functions admitting an infinite amount of roots.","category":"page"},{"location":"internals/","page":"Internals","title":"Internals","text":"Two strategies are currently available: a breadth-first strategy (leaves closer to the root of the tree are processed first); and a depth-first strategy (leaves further away from the root are processed first).","category":"page"},{"location":"internals/#Contractors","page":"Internals","title":"Contractors","text":"","category":"section"},{"location":"internals/","page":"Internals","title":"Internals","text":"To determine the status of a region, the algorithm uses so-called contractors. When we contract a region (wrapped in a Root object), it returns the status of the root and the region. The contractors are various methods to guarantee and refine the status of a root. The available contractors are Bisection, Newton or Krawczyk.","category":"page"},{"location":"internals/","page":"Internals","title":"Internals","text":"julia> using IntervalArithmetic.Symbols\n\njulia> contract(Newton, sin, cos, Root(pi ± 0.001, :unknown))\nRoot([3.14159, 3.1416]_com, :unique)\n\njulia> contract(Newton, sin, cos, Root(2 ± 0.001, :unknown))\nRoot([1.99899, 2.00101]_com, :empty)","category":"page"},{"location":"internals/#RootProblem-and-search-object","page":"Internals","title":"RootProblem and search object","text":"","category":"section"},{"location":"internals/","page":"Internals","title":"Internals","text":"The parameters of a search are represented by a RootProblem object that has the same signature as the roots function.","category":"page"},{"location":"internals/","page":"Internals","title":"Internals","text":"The RootProblem can be iterated to perform the search of the roots, for example to log information at each iteration.","category":"page"},{"location":"internals/","page":"Internals","title":"Internals","text":"For example, the following stops the search after 15 iterations and shows the state of the search at that point.","category":"page"},{"location":"internals/","page":"Internals","title":"Internals","text":"julia> f(x) = exp(x) - sin(x)\nf (generic function with 1 method)\n\njulia> problem = RootProblem(f, interval(-10, 10))\n\njulia> state = nothing   # stores current state of the search\n\njulia> for (k, s) in enumerate(problem)\n           global state = s\n           k == 15 && break  # stop at iteration 15\n       end\n\njulia> state.tree\nBranching\n└─ Branching\n   ├─ Branching\n   │  ├─ Branching\n   │  │  ├─ (:working, Root([-10.0, -8.7886]_com, :unknown))\n   │  │  └─ (:working, Root([-8.78861, -7.55813]_com, :unknown))\n   │  └─ (:final, Root([-6.28132, -6.28131]_com, :unique))\n   └─ Branching\n      ├─ (:working, Root([-5.07783, -3.84734]_com, :unknown))\n      └─ (:working, Root([-3.84736, -2.5975]_com, :unknown))","category":"page"},{"location":"internals/","page":"Internals","title":"Internals","text":"The elements of the iteration are SearchState from the BranchAndPrune.jl package. In the example, we show the tree that get constructed during the search, which, at iteration 15, has found one root and have 4 regions of unknown status to process.","category":"page"}]
}
