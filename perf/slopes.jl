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

struct SlopesMulti
    f::Function
    x::IntervalBox
    c::Vector
    sol::Matrix{Interval}
end

function benchmark_multi()

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
        println("\nFunction $i")
        println("Slope Expansion: ")
        println(DataFrame(slope(rts[i].f, rts[i].x, rts[i].c)))
        println("\nJacobian: ")
        println(DataFrame(ForwardDiff.jacobian(rts[i].f, rts[i].x)))
    end
end
