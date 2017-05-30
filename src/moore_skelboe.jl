# Moore--Skelboe algorithm for global optimization

export moore_skelboe

const inf=infimum
const sup=supremum

"""
    moore_skelboe(X, f, tol=1e-3)

Moore-Skelboe algorithm for global optimization.

"""
function moore_skelboe{T}(f, X::T, tol=1e-3)

    working = [X]
    global_minimum = ∞
    minimizers = []

    while !isempty(working)
        # @show working
        X = pop!(working)
        Y = f(X)

        # @show minimizers

        if inf(Y) < global_minimum

            m = sup(f(Interval(mid(X))))

            if m < global_minimum
                global_minimum = m
            end

            if diam(Y) > tol
                push!(working, bisect(X)...)

            else
                push!(minimizers, X)
            end

        end

    end

    return global_minimum, minimizers
end
