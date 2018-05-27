using IntervalArithmetic, IntervalRootFinding
using ForwardDiff, StaticArrays
using Base.Test

struct Slopes{T}
    f::Function
    x::Interval{T}
    c::T
    sol::Interval{T}
end

@testset "Automatic slope expansion(Float64)" begin

    s = interval(0.75, 1.75)
    rts = Slopes{Float64}[]
    push!(rts, Slopes(x->((x + sin(x)) * exp(-x^2)), s, mid(s), interval(-2.8, 0.1)))
    push!(rts, Slopes(x->(x^4 - 10x^3 + 35x^2 - 50x + 24), s, mid(s), interval(-44, 38.5)))
    push!(rts, Slopes(x->((log(x + 1.25) - 0.84x) ^ 2), s, mid(s), interval(-0.16, 0.44)))
    push!(rts, Slopes(x->(0.02x^2 - 0.03exp(-(20(x - 0.875))^2)), s, mid(s), interval(0.03, 0.33)))
    push!(rts, Slopes(x->(exp(x^2)), s, mid(s), interval(6.03, 33.23)))
    push!(rts, Slopes(x->(x^4 - 12x^3 + 47x^2 - 60x - 20exp(-x)), s, mid(s), interval(-39, 65.56)))
    push!(rts, Slopes(x->(x^6 - 15x^4 + 27x^2 + 250), s, mid(s), interval(-146.9, 67.1)))
    push!(rts, Slopes(x->(atan(cos(tan(x)))), s, mid(s), interval(1, 2)))
    push!(rts, Slopes(x->(asin(cos(acos(sin(x))))), s, mid(s), interval(1.36, ∞)))

    for i in 1:length(rts)
            @test slope(rts[i].f, rts[i].x, rts[i].c) ⊆ rts[i].sol
    end
end

@testset "Automatic slope expansion(BigFloat)" begin
    s = interval(BigFloat(0.75), BigFloat(1.75))
    rts = Slopes{BigFloat}[]
    push!(rts, Slopes(x->((x + sin(x)) * exp(-x^2)), s, mid(s), interval(BigFloat(-2.8), BigFloat(0.1))))
    push!(rts, Slopes(x->(x^4 - 10x^3 + 35x^2 - 50x + 24), s, mid(s), interval(BigFloat(-44), BigFloat(38.5))))
    push!(rts, Slopes(x->((log(x + 1.25) - 0.84x) ^ 2), s, mid(s), interval(BigFloat(-0.16), BigFloat(0.44))))
    push!(rts, Slopes(x->(0.02x^2 - 0.03exp(-(20(x - 0.875))^2)), s, mid(s), interval(BigFloat(0.03), BigFloat(0.33))))
    push!(rts, Slopes(x->(exp(x^2)), s, mid(s), interval(BigFloat(6.03), BigFloat(33.23))))
    push!(rts, Slopes(x->(x^4 - 12x^3 + 47x^2 - 60x - 20exp(-x)), s, mid(s), interval(BigFloat(-39), BigFloat(65.56))))
    push!(rts, Slopes(x->(x^6 - 15x^4 + 27x^2 + 250), s, mid(s), interval(BigFloat(-146.9), BigFloat(67.1))))
    push!(rts, Slopes(x->(atan(cos(tan(x)))), s, mid(s), interval(BigFloat(1), BigFloat(2))))
    push!(rts, Slopes(x->(asin(cos(acos(sin(x))))), s, mid(s), interval(BigFloat(1.36), BigFloat(∞))))

    for i in 1:length(rts)
            @test slope(rts[i].f, rts[i].x, rts[i].c) ⊆ rts[i].sol
    end
end

struct SlopesMulti
    f::Function
    x::IntervalBox
    c::Vector
    sol::Matrix{Interval}
end

@testset "Multidim slope expansion" begin

    rts = SlopesMulti[]
    f(x, y) = SVector(x^2 + y^2 - 1, y - 2x)
    f(X) = f(X...)
    X = (-6..6) × (-6..6)
    c = [0.0, 0.0]
    push!(rts, SlopesMulti(f, X, c, [-6..6 -6..6; -2.. -2 1..1]))

    function g(x)
        (x1, x2, x3) = x
        SVector(    x1^2 + x2^2 + x3^2 - 1,
                    x1^2 + x3^2 - 0.25,
                    x1^2 + x2^2 - 4x3
                )
    end

    X = (-5..5)
    XX = IntervalBox(X, 3)
    cc = [0.0, 0.0, 0.0]
    push!(rts, SlopesMulti(g, XX, cc, [-5..5 -5..5 -5..5; -5..5 0..0 -5..5; -5..5 -5..5 -4.. -4]))
    function h(x)
        (x1, x2, x3) = x
        SVector(    x1 + x2^2 + x3^2 - 1,
                    x1^2 + x3 - 0.25,
                    x1^2 + x2 - 4x3
                )
    end

    XXX = IntervalBox(1..5, 2..6, -3..7)
    ccc = [3.0, 4.0, 2.0]
    push!(rts, SlopesMulti(h, XXX, ccc, [1..1 6..10 -1..9; 4..8 0..0 1..1; 4..8 1..1 -4.. -4]))

    for i in 1:length(rts)
            @test slope(rts[i].f, rts[i].x, rts[i].c) == rts[i].sol
    end

end
