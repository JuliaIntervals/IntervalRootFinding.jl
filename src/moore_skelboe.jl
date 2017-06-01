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

    while !isempty(working)

        X = pop!(working)
        Y = f(X)

        if inf(Y) <= global_minimum

            m = sup(f(Interval.(mid.(X))))

            if m < global_minimum
                global_minimum = m

                # # remove mimimizers that are now excluded:
                #
                # new_minimizers =


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


"""
Do binary search for item `x` in *sorted* vector `v`.
Returns the lower bound for the position of `x` in `v`.
"""
function binary_search(v, x)
    a, b = 1, length(v)
    m = (a + b) ÷ 2  # mid-point

    while abs(a - b) > 1

        m = (a + b) ÷ 2  # mid-point

        if v[m] == x
            return m
        end

        if x < v[m]
            b = m
        else
            a = m
        end
    end

    if v[b] == x
        return b
    end

    return a
end

v = [1, 3, 6, 7, 9, 10]
binary_search(v, 3)
binary_search(v, 4)
