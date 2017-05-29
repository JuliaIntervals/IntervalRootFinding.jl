using IntervalRootFinding
Root = IntervalRootFinding.Root

import IntervalArithmetic.diam

diam(x::Root) = diam(x.interval)

export branch_and_prune, bisection_helper, newton_helper

Base.size(x::Interval) = (1,)

"""
Generic branch and prune.

Inputs:
- Interval or IntervalBox X
- function `prune` that takes as argument an interval and returns
a pair (Symbol, Interval); the Symbol describes the status of that interval
"""
function branch_and_prune(X, f, prune_helper, tol=1e-3)

    prune = prune_helper(f)

    input_dim = size(X)[1]
    output_dim = size(f(X))[1]

    @show input_dim
    @show output_dim

    working = [X]
    outputs = Root{typeof(X)}[]

    while !isempty(working)
        # @show working
        X = pop!(working)

        status, output = prune(X)

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


function bisection_helper(f)
    X -> begin
        image = f(X)

        if !(zero(X) ⊆ image)
            return :empty, X
        end

        return :unknown, X
    end
end

"""
Generic refine operation for Krawczyk and Newton.
This assumes that it is already known that `X` contains a unique root.
Call using e.g. `op = X -> N(f, f_prime, X)`
"""
function refine(op, X::Interval)

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
    f -> begin

        f_prime = x -> ForwardDiff.derivative(f, x)

        X -> begin

            NX = op(f, f_prime, X) ∩ X

            isempty(NX) && return (:empty, X)

            if NX ⪽ X  # isinterior; know there's a unique root inside
                NX =  refine(X -> op(f, f_prime, X), NX)
                return :unique, NX
            end

            return :unknown, NX

        end
    end
end
