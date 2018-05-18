using IntervalArithmetic, IntervalRootFinding

struct Quad{T}
    a::Interval{T}
    b::Interval{T}
    c::Interval{T}
    x::Vector{Interval{T}}
end


@testset "Quadratic Interval Equations" begin
    rts = Quad{Float64}[]
    push!(rts, Quad(interval(1, 5), interval(2, 4), interval(0, 2), [interval(-4, -2), interval(-0, -0)]))
    push!(rts, Quad(interval(-1, 1), interval(-2, 2), interval(-1, 1), [interval(-∞, -1), interval(-1, -1), interval(-1, ∞)]))
    push!(rts, Quad(interval(1, 2), interval(1, 3), interval(2, 6), [interval(-2, -1)]))

    for i in 1:length(rts)
            @test quadratic_roots(rts[i].a, rts[i].b, rts[i].c) == rts[i].x
    end
end
