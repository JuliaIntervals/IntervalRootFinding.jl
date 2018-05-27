# Reference : Dietmar Ratz - An Optimized Interval Slope Arithmetic and its Application
using IntervalArithmetic, ForwardDiff, StaticArrays
import Base: +, -, *, /, ^, sqrt, exp, log, sin, cos, tan, asin, acos, atan
import IntervalArithmetic: mid, interval

"""
Expands the slope of `f` over `x` at the point `c` (default `c = mid(x)`)
"""
function slope(f::Function, x::Interval, c::Number = mid(x))
    try
        f(slope_var(x, c)).s
    catch y
        if isa(y, MethodError)
            ForwardDiff.derivative(f, x)
        end
    end
end

function slope(f::Function, x::Union{IntervalBox, SVector}, c::AbstractVector = mid.(x)) # multidim
    try
        T = typeof(x[1].lo)
        n = length(x)
        s = Vector{Slope{T}}[]
        for i in 1:n
            arr = fill(Slope(zero(T)), n)
            arr[i] = slope_var(x[i], c[i])
            push!(s, arr)
        end
        return slope.(hcat(([(f(s[i])) for i=1:n])...))
    catch y
        if isa(y, MethodError)
            ForwardDiff.jacobian(f, x)
        end
    end

end

struct Slope{T}
    x::Interval{T}
    c::Interval{T}
    s::Interval{T}
end

Slope(c) = Slope(c, c, 0)
Slope(a, b, c) = Slope(promote(convert(Interval, a), b, c)...)

function slope_var(v::Number)
    Slope(v, v, 1)
end

function slope_var(v::Interval, c::Number)
    Slope(v, c, 1)
end

function interval(u::Slope)
    u.x
end

function mid(u::Slope)
    u.c
end

function slope(u::Slope)
    u.s
end

function +(u::Slope, v::Slope)
    Slope(u.x + v.x, u.c + v.c, u.s + v.s)
end

function -(u::Slope, v::Slope)
    Slope(u.x - v.x, u.c - v.c, u.s - v.s)
end

function *(u::Slope, v::Slope)
    Slope(u.x * v.x, u.c * v.c, u.s * v.c + u.x * v.s)
end

function /(u::Slope, v::Slope)
    Slope(u.x / v.x, u.c / v.c, (u.s - (u.c / v.c) * v.s) / v.x)
end

function +(u, v::Slope)
    Slope(u + v.x, u + v.c, v.s)
end

function -(u, v::Slope)
    Slope(u - v.x, u - v.c, -v.s)
end

function *(u, v::Slope)
    Slope(u * v.x, u * v.c, u * v.s)
end

function /(u, v::Slope)
    Slope(u / v.x, u / v.c, -(u / v.c) * (v.s / v.x))
end

+(v::Slope, u) = u + v

*(v::Slope, u) = u * v

function -(u::Slope, v)
    Slope(u.x - v, u.c - v, u.s)
end

function -(u::Slope)
    Slope(-u.x, -u.c, -u.s)
end

function /(u::Slope, v)
    Slope(u.x / v, u.c / v, u.s / v)
end

function sqr(u::Slope)
    Slope(u.x ^ 2, u.c ^ 2, (u.x + u.c) * u.s)
end

function ^(u::Slope, k::Integer)
    if k == 0
        return Slope(1)
    elseif k == 1
        return u
    elseif k == 2
        return sqr(u)
    else
        hxi = interval(u.x.lo) ^ k
        hxs = interval(u.x.hi) ^ k
        hx = hull(hxi, hxs)

        if (k % 2 == 0) && (0 ∈ u.x)
            hx = interval(0, hx.hi)
        end

        hc = u.c ^ k

        i = u.x.lo - u.c.lo
        s = u.x.hi - u.c.hi

        if ((i == 0) || (s == 0) || (k % 2 == 1 && zero(u.x) ⪽ u.x))
            h1 = k * (u.x ^ (k - 1))
        else
            if k % 2 == 0 || u.x.lo >= 0
                h1 = interval((hxi.hi - hc.lo) / i, (hxs.hi - hc.lo) / s)
            else
                h1 = interval((hxs.lo - hc.hi) / s, (hxi.lo - hc.hi) / i)
            end
        end
        return Slope(hx, hc, h1 * u.s)
    end
end

function sqrt(u::Slope)
    Slope(sqrt(u.x), sqrt(u.c), u.s / (sqrt(u.x) + sqrt(u.c)))
end

function exp(u::Slope)
    hx = exp(u.x)
    hc = exp(u.c)

    i = u.x.lo - u.c.lo
    s = u.x.hi - u.c.hi

    if (i == 0 || s == 0)
        h1 = hx
    else
        h1 = interval((hx.lo - hc.lo) / i, (hx.hi - hc.hi) / s)
    end

    Slope(hx, hc, h1 * u.s)
end

function log(u::Slope)
    hx = log(u.x)
    hc = log(u.c)

    i = u.x.lo - u.c.lo
    s = u.x.hi - u.c.hi

    if (i == 0 || s == 0)
        h1 = 1 / u.x
    else
        h1 = interval((hx.hi - hc.hi) / s, (hx.lo - hc.lo) / i)
    end
    Slope(hx, hc, h1 * u.s)
end

function sin(u::Slope) # Using derivative to upper bound the slope expansion for now
    hx = sin(u.x)
    hc = sin(u.c)
    hs = cos(u.x)
    Slope(hx, hc, hs)
end

function cos(u::Slope) # Using derivative to upper bound the slope expansion for now
    hx = cos(u.x)
    hc = cos(u.c)
    hs = -sin(u.x)
    Slope(hx, hc, hs)
end

function tan(u::Slope) # Using derivative to upper bound the slope expansion for now
    hx = tan(u.x)
    hc = tan(u.c)
    hs = (sec(u.x)) ^ 2
    Slope(hx, hc, hs)
end

function asin(u::Slope)
    hx = asin(u.x)
    hc = asin(u.c)
    hs = 1 / sqrt(1 - (u.x ^ 2))
    Slope(hx, hc, hs)
end

function acos(u::Slope)
    hx = acos(u.x)
    hc = acos(u.c)
    hs = -1 / sqrt(1 - (u.x ^ 2))
    Slope(hx, hc, hs)
end

function atan(u::Slope)
    hx = atan(u.x)
    hc = atan(u.c)
    hs = 1 / 1 + (u.x ^ 2)
    Slope(hx, hc, hs)
end
