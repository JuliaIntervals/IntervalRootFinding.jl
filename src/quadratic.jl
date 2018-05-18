"""
Helper function for quadratic_interval
"""
function quadratic_helper!{T}(a::Interval{T}, b::Interval{T}, c::Interval{T}, L::Array{Interval{T}})

    Δ = b*b - 4*a*c

    Δ.hi < 0 && return

    Δ = sqrt(Δ)

    if (b.lo >= 0)
        x0 = -0.5 * (b + Δ)

    else
        x0 = -0.5 * (b - Δ)
    end

    if (0 ∈ x0)
        push!(L, x0)

    else
        x1 = c / x0
        x0 = x0 / a
        push!(L, x0, x1)
    end

end


"""
Function to solve a quadratic equation where the coefficients are intervals.
Returns an array of intervals of the roots.
Eldon Hansen and G. William Walster : Global Optimization Using Interval Analysis - Chapter 8
"""
function quadratic_roots{T}(a::Interval{T}, b::Interval{T}, c::Interval{T})

    L = Interval{T}[]
    R = Interval{T}[]

    quadratic_helper!(Interval(a.lo), Interval(b.lo), Interval(c.lo), L)
    quadratic_helper!(Interval(a.hi), Interval(b.hi), Interval(c.hi), L)
    quadratic_helper!(Interval(a.lo), Interval(b.hi), Interval(c.lo), L)
    quadratic_helper!(Interval(a.hi), Interval(b.lo), Interval(c.hi), L)

    if (length(L) == 8)
        resize!(L, 4)
    end

    if (a.lo < 0 || (a.lo == 0 && b.hi == 0) || (a.lo == 0 && b.hi == 0 && c.lo ≤ 0))
        push!(L, Interval(-∞))
    end

    if (a.lo < 0 || (a.lo == 0 && b.lo == 0) || (a.lo == 0 && b.lo == 0 && c.lo ≤ 0))
        push!(L, Interval(∞))
    end

    sort!(L, by = x -> x.lo)

    for i in 1:2:length(L)
        push!(R, Interval(L[i].lo, L[i+1].hi))
    end

    R
end
