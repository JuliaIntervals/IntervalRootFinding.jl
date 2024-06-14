module IntervalRootFinding

using Reexport
@reexport using IntervalArithmetic
using IntervalArithmetic.Symbols

using BranchAndPrune
using ForwardDiff
using ForwardDiff: derivative, jacobian
using StaticArrays

using LinearAlgebra

import Base: ⊆, show, big, \

export derivative, jacobian

include("region.jl")
export isempty_region, intersect_region, in_region

include("root_object.jl")
export Root, isunique, root_status, root_region

include("roots.jl")
export roots, RootProblem

include("contractors.jl")
export Bisection, Newton, Krawczyk, contract

include("linear_eq.jl")
export gauss_seidel_interval, gauss_elimination_interval, gauss_seidel_contractor

include("complex.jl")

include("newton1d.jl")
include("quadratic.jl")
include("slopes.jl")

end