# Moore--Skelboe algorithm for global optimization

export moore_skelboe

const inf=infimum
const sup=supremum

"""
    moore_skelboe(f, X, tol=1e-3)

Moore-Skelboe algorithm for global optimization.

"""
function moore_skelboe{T}(f, X::T, tol=1e-3)

    working = [X]
    global_minimum = ∞
    minimizers = T[]

    num_bisections = 0

    while !isempty(working)

        X = pop!(working)
        Y = f(X)

        if inf(Y) > global_minimum
            continue
        end

        m = sup(f(Interval.(mid.(X))))

        if m < global_minimum
            global_minimum = m
        end

        if diam(Y) > tol
            push!(working, bisect(X)...)
            num_bisections += 1

        else
            push!(minimizers, X)
        end

    end

    @show num_bisections

    return global_minimum, minimizers
end


using DataStructures
export optimize

"""
    optimize(f, X, tol=1e-3)

"Optimize" algorithm for global optimization, not including constraint propagation.
Removes boxes that are known not to contain the global minimum

"""
function optimize{T}(f, X::T, tol=1e-3)

    # list of boxes with corresponding lower bound, ordered by increasing lower bound:
    working = SortedVector([(X, ∞)], x->x[2])

    global_min = ∞
    minimizers = Tuple{T, Float64}[]

    num_bisections = 0

    while !isempty(working)

        # @show working

        (X, Xmin) = pop!(working)
        # Y = f(X)

        if inf(f(X)) > global_min  # inf(f(X)) is Xmin?
            continue
        end

        # find candidate for global minimum upper bound:
        m = sup(f(Interval.(mid.(X))))   # evaluate at midpoint of current interval

        # @show m, global_min
        # @show length(working.data)

        if m < global_min
            global_min = m
        end

        # Remove all boxes whose lower bound is greater than the current one:
        # Since they are ordered, just find the first one that is too big

        cutoff = searchsortedfirst(working.data, (X, global_min), by=x->x[2])
        resize!(working, cutoff-1)




        if diam(X) < tol
            push!( minimizers, (X, inf(f(X))) )

        else
            X1, X2 = bisect(X)
            push!( working, (X1, inf(f(X1))), (X2, inf(f(X2))) )
            num_bisections += 1
        end

    end

    return global_min, minimizers
end
