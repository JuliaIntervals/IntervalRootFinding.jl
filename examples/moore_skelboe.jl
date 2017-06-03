using IntervalArithmetic, IntervalRootFinding

moore_skelboe(x->(x-1)^2*(x+1)^2 + 1, -5..5)
optimize(x->(x-1)^2*(x+1)^2 + 1, -5..5)

# Branin function (Jaulin pg. 120)

a = @interval(5.1 / (4*(pi^2)))
b = @interval(5 / pi)
c = @interval(1 - (1 / (8pi)))

f(p) = ( (p1, p2) = p; pow(p2 - a*pow(p1,2) + b*p1 - 6, 2) +10*c * cos(p1) +10 )

@time global_min, minimizers = moore_skelboe(f, (-5..10)×(0..15))

@time global_min, minimizers = optimize(f, (-5..10)×(0..15))


minimizers
