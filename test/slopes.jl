using IntervalArithmetic, IntervalRootFinding
using ForwardDiff
using Test

struct Slopes{T}
    f::Function
    x::Interval{T}
    c::T
    sol::Interval{T}
end

@testset "Automatic slope expansion" begin
    for T in [Float64, BigFloat]
        s = interval(T(0.75), T(1.75))
        example = Slopes{T}[]
        push!(example, Slopes(x->((x + sin(x)) * exp(-x^2)), s, mid(s), interval(T(-2.8), T(0.1))))
        push!(example, Slopes(x->(x^4 - 10x^3 + 35x^2 - 50x + 24), s, mid(s), interval(T(-44), T(38.5))))
        push!(example, Slopes(x->((log(x + 1.25) - 0.84x) ^ 2), s, mid(s), interval(T(-0.16), T(0.44))))
        push!(example, Slopes(x->(0.02x^2 - 0.03exp(-(20(x - 0.875))^2)), s, mid(s), interval(T(0.03), T(0.33))))
        push!(example, Slopes(x->(exp(x^2)), s, mid(s), interval(T(6.03), T(33.23))))
        push!(example, Slopes(x->(x^4 - 12x^3 + 47x^2 - 60x - 20exp(-x)), s, mid(s), interval(T(-39), T(65.56))))
        push!(example, Slopes(x->(x^6 - 15x^4 + 27x^2 + 250), s, mid(s), interval(T(-146.9), T(67.1))))
        push!(example, Slopes(x->(atan(cos(tan(x)))), s, mid(s), interval(T(1), T(2))))
        push!(example, Slopes(x->(asin(cos(acos(sin(x))))), s, mid(s), interval(T(1.36), T(Inf))))

        for i in eachindex(example)
            s = slope(example[i].f, example[i].x, example[i].c)
            @test issubset_interval(s, example[i].sol)
        end
    end
end
