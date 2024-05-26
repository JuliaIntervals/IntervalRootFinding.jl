# This file is part of the ValidatedNumerics.jl package; MIT licensed

module IntervalRootFinding

using Reexport
@reexport using IntervalArithmetic
using IntervalArithmetic.Symbols

using BranchAndPrune
using ForwardDiff
using StaticArrays

using ForwardDiff: derivative, gradient, jacobian
using LinearAlgebra

import Base: ⊆, show, big, \

## Root finding
export
    derivative, jacobian,  # reexport derivative from ForwardDiff
    Root, isunique, root_status, root_region,
    roots, newton1d, quadratic_roots,
    gauss_seidel_interval, gauss_seidel_interval!,
    gauss_seidel_contractor, gauss_seidel_contractor!,
    gauss_elimination_interval, gauss_elimination_interval!,
    slope

import IntervalArithmetic: interval

include("region.jl")
export
    isempty_region, intersect_region, in_region

include("root_object.jl")
include("roots.jl")

include("complex.jl")
include("contractors.jl")
include("newton1d.jl")
include("quadratic.jl")
include("linear_eq.jl")
include("slopes.jl")

end