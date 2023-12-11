# This file is part of the ValidatedNumerics.jl package; MIT licensed

module IntervalRootFinding

using Reexport
@reexport using IntervalArithmetic
using IntervalArithmetic.Symbols

using BranchAndPrune
using ForwardDiff
using StaticArrays

using ForwardDiff: derivative, gradient, jacobian
using LinearAlgebra: I, Diagonal

import Base: ⊆, show, big, \

## Root finding
export
    derivative, jacobian,  # reexport derivative from ForwardDiff
    Root, isunique, root_status,
    roots, find_roots,
    bisect, newton1d, quadratic_roots,
    gauss_seidel_interval, gauss_seidel_interval!,
    gauss_seidel_contractor, gauss_seidel_contractor!,
    gauss_elimination_interval, gauss_elimination_interval!,
    slope

import IntervalArithmetic: interval


const D = derivative


const where_bisect = 0.49609375  # 127//256

const Region = Union{Interval, SVector{N, Interval} where N}

IntervalBox(x::Interval, N::Integer) = SVector{N}(fill(x, N))
IntervalBox(xx::Vararg{Interval, N}) where N = SVector{N}(xx...)

include("root_object.jl")

include("complex.jl")
include("contractors.jl")
include("roots.jl")
include("newton1d.jl")
include("quadratic.jl")
include("linear_eq.jl")
include("slopes.jl")

end