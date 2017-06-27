using IntervalRootFinding
Root = IntervalRootFinding.Root

import IntervalArithmetic: diam, isinterior

diam(x::Root) = diam(x.interval)

export branch_and_prune, bisection_helper, newton_helper

Base.size(x::Interval) = (1,)

isinterior{N}(X::IntervalBox{N}, Y::IntervalBox{N}) = all(isinterior.(X, Y))


"""
    branch_and_prune(X, f, contractor, tol=1e-3)

Generic branch and prune routine for finding isolated roots of a function ``f:R^n → R^n`` in a box.

Inputs:
- `X`: `Interval` or `IntervalBox`
- `f`: function whose roots will be found
- `contractor`: function that, when applied to the function `f`, determines the status of a given box `X`. It returns the new box and a symbol indicating the status.
"""
function branch_and_prune(X, f, contractor, tol=1e-3)

    input_dim = size(X)[1]
    output_dim = size(f(X))[1]

    # @show input_dim
    # @show output_dim

    # if !(input_dim == output_dim)
    #     throw(ArgumentError("Input dimension ($input_dim) and output dimension ($output_dim) must be the same."))
    # end

    contract = contractor(input_dim, f)

    # main algorithm:

    working = [X]
    outputs = Root{typeof(X)}[]

    sizehint!(outputs, 100)
    sizehint!(working, 1000)

    while !isempty(working)
        # @show working
        X = pop!(working)

        status, output = contract(X)

        if status == :empty
            continue

        elseif status == :unique
            push!(outputs, Root(output, :unique))

        elseif diam(output) < tol
            push!(outputs, Root(output, :unknown))

        else  # branch
            X1, X2 = bisect(X)

            push!(working, X1, X2)
        end

    end

    return outputs
end


abstract type Contractor end

struct BisectionContractor{F} <: Contractor
    dimension::Int
    f::F
end

contains_zero{T}(X::Interval{T}) = zero(T) ∈ X
contains_zero(X::IntervalBox) = all(contains_zero(X[i]) for i in 1:length(X))

function (contractor::BisectionContractor)(X)
    image = contractor.f(X)

    if !(contains_zero(image))
        return :empty, X
    end

    return :unknown, X
end



"""
Generic refine operation for Krawczyk and Newton.
This assumes that it is already known that `X` contains a unique root.
Call using e.g. `op = X -> N(f, f_prime, X)`
"""
function refine(op, X)

    tolerance = 1e-16

    while diam(X) > tolerance  # avoid problem with tiny floating-point numbers if 0 is a root
        NX = op(X) ∩ X
        NX == X && break  # reached limit of precision
        X = NX
    end

    return X
end

"""
Helper function for Newton (`op=N`) and Krawczyk (`op=K`)
"""
function newton_helper(op)
    (dim, f) -> begin

        if dim == 1
            f_prime = x -> ForwardDiff.derivative(f, x)
        else
            f_prime = x -> ForwardDiff.jacobian(f, x)
        end

        X -> begin

            if !(zero(X) ⊆ IntervalBox(f(X)...))
                return :empty, X
            end


            NX = op(f, f_prime, X) ∩ X

            isempty(NX) && return :empty, X


            if NX ⪽ X  # isinterior; know there's a unique root inside
                NX =  refine(X -> op(f, f_prime, X), NX)
                return :unique, NX
            end


            return :unknown, NX

        end
    end
end
