using BenchmarkTools, Compat, DataFrames, IntervalRootFinding, IntervalArithmetic, StaticArrays

function benchmark_results()
    f = [] # Example functions from Dietmar Ratz - An Optimized Interval Slope Arithmetic and its Application
    push!(f, x->((x + sin(x)) * exp(-x^2)))
    push!(f, x->(x^4 - 10x^3 + 35x^2 - 50x + 24))
    push!(f, x->((log(x + 1.25) - 0.84x) ^ 2))
    push!(f, x->(0.02x^2 - 0.03exp(-(20(x - 0.875))^2)))
    push!(f, x->(exp(x^2)))
    push!(f, x->(x^4 - 12x^3 + 47x^2 - 60x - 20exp(-x)))
    push!(f, x->(x^6 - 15x^4 + 27x^2 + 250))
    push!(f, x->(atan(cos(tan(x)))))
    push!(f, x->(asin(cos(acos(sin(x))))))

    s = interval(0.75, 1.75)
    df = DataFrame()

    df[:Method] = ["Automatic Differentiation", "Slope Expansion"]
    for n in 1:length(f)

        t1 = ForwardDiff.derivative(f[n], s)
        t2 = slope(f[n], s, mid(s))
        df[Symbol("f" * "$n")] = [t1, t2]
    end
    a = []
    for i in 1:length(f)
        push!(a, Symbol("f" * "$i"))
    end
    df1 = stack(df, a)
    dfnew = unstack(df1, :variable, :Method, :value)
    dfnew = rename(dfnew, :variable => :Function)
    println(dfnew)
    dfnew
end

function benchmark_time()
    f = []
    push!(f, x->((x + sin(x)) * exp(-x^2)))
    push!(f, x->(x^4 - 10x^3 + 35x^2 - 50x + 24))
    push!(f, x->((log(x + 1.25) - 0.84x) ^ 2))
    push!(f, x->(0.02x^2 - 0.03exp(-(20(x - 0.875))^2)))
    push!(f, x->(exp(x^2)))
    push!(f, x->(x^4 - 12x^3 + 47x^2 - 60x - 20exp(-x)))
    push!(f, x->(x^6 - 15x^4 + 27x^2 + 250))
    push!(f, x->(atan(cos(tan(x)))))
    push!(f, x->(asin(cos(acos(sin(x))))))

    s = interval(0.75, 1.75)
    df = DataFrame()
    df[:Method] = ["Automatic Differentiation", "Slope Expansion"]
    for n in 1:length(f)

        t1 = @belapsed ForwardDiff.derivative($f[$n], $s)
        t2 = @belapsed slope($f[$n], $s, mid($s))
        df[Symbol("f" * "$n")] = [t1, t2]
    end
    a = []
    for i in 1:length(f)
        push!(a, Symbol("f" * "$i"))
    end
    df1 = stack(df, a)
    dfnew = unstack(df1, :variable, :Method, :value)
    dfnew = rename(dfnew, :variable => :Function)
    println(dfnew)
    dfnew
end
