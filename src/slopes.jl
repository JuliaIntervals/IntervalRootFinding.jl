# Reference : Dietmar Ratz - An Optimized Interval Slope Arithmetic and its Application
using IntervalArithmetic
import Base: +, -, *, /, ^, sqrt, exp, log

function slope(f::Function, x::Interval, c::Real)
    f(SlopeVar(x, c)).fs
end

struct SlopeType
    fx::Interval
    fc::Interval
    fs::Interval
end

function SlopeConst(c::Union{Real, Interval})
    SlopeType(c, c, 0)
end

function SlopeVar(v::Real)
    SlopeType(v, v, 1)
end

function SlopeVar(v::Interval, c::Real)
    SlopeType(v, c, 1)
end

function fxValue(u::SlopeType)
    u.fx
end

function fcValue(u::SlopeType)
    u.fc
end

function fsValue(u::SlopeType)
    u.fs
end

function +(u::SlopeType, v::SlopeType)
    SlopeType(u.fx + v.fx, u.fc + v.fc, u.fs + v.fs)
end

function -(u::SlopeType, v::SlopeType)
    SlopeType(u.fx - v.fx, u.fc - v.fc, u.fs - v.fs)
end

function *(u::SlopeType, v::SlopeType)
    SlopeType(u.fx * v.fx, u.fc * v.fc, u.fs * v.fc + u.fx * v.fs)
end

function /(u::SlopeType, v::SlopeType)
    SlopeType(u.fx / v.fx, u.fc / v.fc, (u.fs - (u.fc / v.fc) * v.fs) / v.fx)
end

function +(u::Union{Interval, Real}, v::SlopeType)
    SlopeType(u + v.fx, u + v.fc, v.fs)
end

function -(u::Union{Interval, Real}, v::SlopeType)
    SlopeType(u - v.fx, u - v.fc, -v.fs)
end

function *(u::Union{Interval, Real}, v::SlopeType)
    SlopeType(u * v.fx, u * v.fc, u * v.fs)
end

function /(u::Union{Interval, Real}, v::SlopeType)
    SlopeType(u / v.fx, u / v.fc, -(u / v.fc) * (v.fs / v.fx))
end

+(v::SlopeType, u::Union{Interval, Real}) = u + v

-(v::SlopeType, u::Union{Interval, Real}) = u - v

*(v::SlopeType, u::Union{Interval, Real}) = u * v

/(v::SlopeType, u::Union{Interval, Real}) = u / v

function sqr(u::SlopeType)
    SlopeType(u.fx ^ 2, u.fc ^ 2, (u.fx + u.fc) * u.fs)
end

function ^(u::SlopeType, k::Integer)
    if k == 0
        return SlopeConst(1)
    elseif k == 1
        return u
    elseif k == 2
        return sqr(u)
    else
        hxi = interval(u.fx.lo) ^ k
        hxs = interval(u.fx.hi) ^ k
        hx = hull(hxi, hxs)

        if (k % 2 == 0) && (0 ∈ u.fx)
            hx = interval(0, hx.hi)
        end

        hc = u.fc ^ k

        i = u.fx.lo - u.fc.lo
        s = u.fx.hi - u.fc.hi

        if ((i == 0) || (s == 0) || (k % 2 == 1 && Interval(0) ⪽ u.fx))
            h1 = k * (u.fx ^ (k - 1))
        else
            if k % 2 == 0 || u.fx.lo >= 0
                h1 = interval((hxi.hi - hc.lo) / i, (hxs.hi - hc.lo) / s)
            else
                h1 = interval((hxs.lo - hc.hi) / s, (hxi.lo - hc.hi) / i)
            end
        end
        return SlopeType(hx, hc, h1 * u.fs)
    end
end

function sqrt(u::SlopeType)
    SlopeType(sqrt(u.fx), sqrt(u.fc), u.fs / (sqrt(u.fx) + sqrt(u.fc)))
end

function exp(u::SlopeType)
    hx = exp(u.fx)
    hc = exp(u.fc)

    i = u.fx.lo - u.fc.lo
    s = u.fx.hi - u.fc.hi

    if (i == 0 || s == 0)
        h1 = hx
    else
        h1 = interval((hx.lo - hc.lo) / i, (hx.hi - hc.hi) / s)
    end

    SlopeType(hx, hc, h1 * u.fs)
end

function log(u::SlopeType)
    hx = log(u.fx)
    hc = log(u.fc)

    i = u.fx.lo - u.fc.lo
    s = u.fx.hi - u.fc.hi

    if (i == 0 || s == 0)
        h1 = 1 / u.fx
    else
        h1 = interval((hx.hi - hc.hi) / s, (hx.lo - hc.lo) / i)
    end
    SlopeType(hx, hc, h1 * u.fs)
end
