using IntervalRootFinding, IntervalArithmetic


"""
Helper function for quadratic_interval
"""
function quadratic_helper(a::Number, b::Number, c::Number, L::AbstractArray)

    Δ = b*b - 4*a*c

    Δ < 0 && return

    Δ = sqrt(Δ)

    if (b >= 0)
        x0 = -0.5 * (b + Δ)

    else
        x0 = -0.5 * (b - Δ)
    end

    if (x0 == 0)
        x1 = x0

    else
        x1 = c / x0
        x0 = x0 / a
        push!(L, x1)
    end

    push!(L, x0)
end


"""
Function to solve a quadratic equation where the coefficients are intervals.
Returns an array of intervals of the roots.
Eldon Hansen and G. William Walster : Global Optimization Using Interval Analysis - Chapter 8
"""
function quadratic_interval{T}(a::Interval{T}, b::Interval{T}, c::Interval{T})

    L = T[]
    R = Interval{T}[]

    quadratic_helper(a.lo, b.lo, c.lo, L)
    quadratic_helper(a.hi, b.hi, c.hi, L)
    quadratic_helper(a.lo, b.hi, c.lo, L)
    quadratic_helper(a.hi, b.lo, c.hi, L)

    if (length(L) == 8)
        resize!(L, 4)
    end

    if (a.lo < 0 || (a.lo == 0 && b.hi == 0) || (a.lo == 0 && b.hi == 0 && c.lo ≤ 0))
        push!(L, -∞)
    end

    if (a.lo < 0 || (a.lo == 0 && b.lo == 0) || (a.lo == 0 && b.lo == 0 && c.lo ≤ 0))
        push!(L, ∞)
    end

    sort!(L)

    for i in 1:2:length(L)
        push!(R, Interval(L[i], L[i+1]))
    end

    R
end
