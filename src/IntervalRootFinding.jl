# This file is part of the ValidatedNumerics.jl package; MIT licensed

__precompile__(true)

module IntervalRootFinding

using IntervalArithmetic
using ForwardDiff
using StaticArrays

using LinearAlgebra: I 


import Base: ⊆, show, big, \

import Polynomials: roots

## Root finding
export
    derivative, jacobian,  # reexport derivative from ForwardDiff
    Root, is_unique,
    roots, find_roots,
    bisect, newton1d, quadratic_roots,
    gauss_seidel_interval, gauss_seidel_interval!,
    gauss_seidel_contractor, gauss_seidel_contractor!,
    gauss_elimination_interval, gauss_elimination_interval!,
    slope

export isunique, root_status


import IntervalArithmetic: interval, wideinterval



const derivative = ForwardDiff.derivative
const D = derivative

const where_bisect = IntervalArithmetic.where_bisect ## 127//256


include("root_object.jl")

include("newton.jl")
include("krawczyk.jl")

include("complex.jl")
include("contractors.jl")
include("roots.jl")
include("newton1d.jl")
include("quadratic.jl")
include("linear_eq.jl")
include("slopes.jl")


gradient(f) = X -> ForwardDiff.gradient(f, X[:])

ForwardDiff.jacobian(f, X::IntervalBox) = ForwardDiff.jacobian(f, X.v)

const ∇ = gradient

export ∇




function find_roots(f::Function, a::Interval{T}, method::Function = newton;
                    tolerance = eps(T), debug = false, maxlevel = 30) where {T}

    method(f, a; tolerance=tolerance, debug=debug, maxlevel=maxlevel)
end

function find_roots(f::Function, f_prime::Function, a::Interval{T}, method::Function=newton;
                    tolerance=eps(T), debug=false, maxlevel=30) where {T}

    method(f, f_prime, a; tolerance=tolerance, debug=debug, maxlevel=maxlevel)
end

function find_roots(f::Function, a::Real, b::Real, method::Function=newton;
           tolerance=eps(1.0*a), debug=false, maxlevel=30, precision::Int=-1)

    if precision >= 0
        setprecision(Interval, precision) do
            find_roots(f, @interval(a, b), method; tolerance=tolerance, debug=debug, maxlevel=maxlevel)
        end


    else  # use current precision

        find_roots(f, @interval(a, b), method; tolerance=tolerance, debug=debug, maxlevel=maxlevel)

    end
end



function find_roots_midpoint(f::Function, a::Real, b::Real, method::Function=newton;
           tolerance=eps(1.0*a), debug=false, maxlevel=30, precision=-1)

    roots = find_roots(f, a, b, method; tolerance=tolerance, debug=debug, maxlevel=maxlevel, precision=precision)

    T = eltype(roots[1].interval)

    midpoints = T[]
    radii = T[]
    root_symbols = Symbol[]  # :unique or :unknown

    if length(roots) == 0
        return (midpoints, radii, root_symbols)  # still empty
    end

    for root in roots
        midpoint, radius = midpoint_radius(root.interval)
        push!(midpoints, midpoint)
        push!(radii, radius)

        push!(root_symbols, root.status)

    end

    (midpoints, radii, root_symbols)

end

function Base.lexcmp(a::Interval{T}, b::Interval{T}) where {T}
    #@show a, b
    if a.lo < b.lo
        return -1
    elseif a.lo == b.lo
        if a.hi < b.hi
            return -1
        else
            return 0
        end
    end
    return 1

end

Base.lexcmp(a::Root{T}, b::Root{T}) where {T} = lexcmp(a.interval, b.interval)


function clean_roots(f, roots)

    # order, remove duplicates, and include intervals X only if f(X) contains 0
    sort!(roots, lt=lexless)
    roots = unique(roots)
    roots = filter(x -> 0 ∈ f(x.interval), roots)

    # merge neighbouring roots if they touch:

    if length(roots) < 2
        return roots
    end


    new_roots = eltype(roots)[]

    base_root = roots[1]

    for i in 2:length(roots)
        current_root = roots[i]

        if isempty(base_root.interval ∩ current_root.interval) ||
                (base_root.status != current_root.status)

            push!(new_roots, base_root)
            base_root = current_root
        else
            base_root = Root(hull(base_root.interval, current_root.interval), base_root.status)
        end
    end

    push!(new_roots, base_root)

    return new_roots

end

end
