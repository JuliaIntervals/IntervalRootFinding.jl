using IntervalArithmetic, IntervalRootFinding

moore_skelboe(x->(x-1)^2*(x+1)^2 + 1, -5..5)
optimize(x->(x-1)^2*(x+1)^2 + 1, -5..5)

# Branin function (Jaulin pg. 120)

a = @interval(5.1 / (4*(pi^2)))
b = @interval(5 / pi)
c = @interval(1 - (1 / (8pi)))

f(p) = ( (p1, p2) = p; pow(p2 - a*pow(p1,2) + b*p1 - 6, 2) +10*c * cos(p1) +10 )

#@time global_min, minimizers = moore_skelboe(f, (-5..10)×(0..15))

@time global_min, minimizers = optimize(f, (-5..10)×(0..15))


minimizers

minimizer_sets = first.(minimizers)

using LightGraphs

n = length(minimizer_sets)
g = Graph()

add_vertices!(g, length(minimizer_sets))

g

for i in 1:n
    for j in i+1:n
        if !isempty(minimizer_sets[i] ∩ minimizer_sets[j])
            add_edge!(g, i, j)
        end
    end
end

g

connected_components(g)


# SIAM example:
f(X) = ((x, y) = X; exp(sin(50x)) + sin(60*exp(y)) + sin(70*sin(x)) + sin(sin(80y)) - sin(10*(x + y)) + (x^2 + y^2) / 4)

@time global_min, minimizers = optimize(f, (-10..10)×(-10..10), 1e-3)

length(minimizers)

images = f.(first.(minimizers))
minimizers
reduce(∪, images)

# Griewank

G(x) = (n = length(x); 1 + sum(pow(x[i],2) for i in 1:n) / 4000 - prod(cos(x[i] / √i) for i in 1:n))

G̃(x) = G(x) + pow(sin(x[1]), 2) + pow(cos(x[1]), 2) - 1

@time optimize(G, IntervalBox(-600..600))

@time optimize(G, (-600..600) × (-600..600))



@time optimize(G, (-600..600) × (-600..600) × (-600..600))



X = -600..600


ntuple

Y(n) = IntervalBox(ntuple(i->-600..600, n))

@time optimize(G, Y(10))

@time optimize(G, Y(40))

@time optimize(G, Y(50))


@time optimize(G̃, Y(1))

@time optimize(G̃, Y(2))

@time optimize(G̃, Y(3))



# Beale function (Schwefel); Hansen & Walster pg. 334
f(x, y) = pow(1.5 − x*(1 − y), 2) + pow(2.25 - x*(1 − pow(y,2)), 2) +
        (2.625 −  x*(1 − pow(y, 3)))^2
