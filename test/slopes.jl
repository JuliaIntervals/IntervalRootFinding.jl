using IntervalArithmetic, IntervalRootFinding
using ForwardDiff
using Base.Test

struct Slopes{T}
    f::Function
    x::Interval{T}
    c::T
    sol::Interval{T}
end

@testset "Automatic slope expansion" begin
    rts = Slopes{Float64}[]
    push!(rts, Slopes(x->((x + sin(x)) * exp(-x^2)), interval(0.75, 1.75), 1.25, interval(-2.8, 0.1)))
    push!(rts, Slopes(x->(x^4 - 10x^3 + 35x^2 - 50x + 24), interval(0.75, 1.75), 1.25, interval(-44, 38.5)))
    push!(rts, Slopes(x->((log(x + 1.25) - 0.84x) ^ 2), interval(0.75, 1.75), 1.25, interval(-0.16, 0.44)))
    push!(rts, Slopes(x->(0.02x^2 - 0.03exp(-(20(x - 0.875))^2)), interval(0.75, 1.75), 1.25, interval(0.04, 0.33)))
    push!(rts, Slopes(x->(exp(x^2)), interval(0.75, 1.75), 1.25, interval(6.03, 33.23)))
    push!(rts, Slopes(x->(x^4 - 12x^3 + 47x^2 - 60x - 20exp(-x)), interval(0.75, 1.75), 1.25, interval(-39, 65.56)))
    push!(rts, Slopes(x->(x^6 - 15x^4 + 27x^2 + 250), interval(0.75, 1.75), 1.25, interval(-146.9, 67.1)))

    for i in 1:length(rts)
            @test slope(rts[i].f, rts[i].x, rts[i].c) ⊆ rts[i].sol
    end
end
