"""`newton1d` performs the interval Newton method on the given function `f`
with its derivative `f′` and initial interval `x`.
Optional keyword arguments give the tolerances `reltol` and `abstol`.
`reltol` is the tolerance on the relative error whereas `abstol` is the tolerance on |f(X)|,
and a `debug` boolean argument that prints out diagnostic information."""

function newton1d{T}(f::Function, f′::Function, x::Interval{T};
                    reltol=eps(T), abstol=eps(T), debug=false, debugroot=false)

    L = Interval{T}[]

    R = Root{Interval{T}}[]
    reps = reps1 = 0

    push!(L, x)
    initial_width = ∞
    X = emptyinterval(T)
    while !isempty(L)
        X = pop!(L)

        debug && (print("Current interval popped: "); @show X)

        m = mid(X)
        if (isempty(X))
            continue
        end

        if 0 ∉ f′(X)

            debug && println("0 ∉ f′(X)")

            while true

                m = mid(X)
                N = m - (f(interval(m)) / f′(X))

                debug && (print("Newton step1: "); @show (X, X ∩ N))
                if X == X ∩ N
                    reps1 += 1
                    if reps1 > 20
                        reps1 = 0
                        break
                    end
                end
                X = X ∩ N

                if (isempty(X) || diam(X) == 0)
                    break

                elseif 0 ∈ f(interval(prevfloat(mid(X)), nextfloat(mid(X))))
                    n = fa = fb = 0
                    root_exist = true
                    while (n < 4 && (fa == 0 || fb == 0))
                        if fa == 0
                            if 0 ∈ f(interval(prevfloat(X.lo), nextfloat(X.lo)))
                                fa = 1
                            else
                                N = X.lo - (f(interval(X.lo)) / f′(X))
                                X = X ∩ N
                                if (isempty(X) || diam(X) == 0)
                                    root_exist = false
                                    break
                                end
                            end
                        end
                        if fb == 0
                            if 0 ∈ f(interval(prevfloat(X.hi), nextfloat(X.hi)))
                                fb = 1
                            else
                                if 0 ∈ f(interval(prevfloat(mid(X)), nextfloat(mid(X))))
                                    N = X.hi - (f(interval(X.hi)) / f′(X))
                                else
                                    N = mid(X) - (f(interval(mid(X))) / f′(X))
                                end
                                X = X ∩ N
                                if (isempty(X) || diam(X) == 0)
                                    root_exist = false
                                    break
                                end
                            end
                        end
                        N = mid(X) - (f(interval(mid(X))) / f′(X))
                        X = X ∩ N
                        if (isempty(X) || diam(X) == 0)
                            root_exist = false
                            break
                        end
                        n += 1
                    end
                    if root_exist
                        push!(R, Root(X, :unique))
                        debugroot && @show "Root found", X
                    end

                    break
                end
            end

        else
            if diam(X) == initial_width
                reps += 1
                if reps > 10
                    push!(R, Root(X, :unknown))
                    debugroot && @show "Repititive root found", X
                    reps = 0
                    continue
                end
            end
            initial_width = diam(X)
            debug && println("0 ∈ f'(X)")

            expansion_pt = Inf
            # expansion point for the newton step might be m, X.lo or X.hi according to some conditions

            if 0 ∈ f(interval(prevfloat(mid(X)), nextfloat(mid(X))))
                # 0 ∈ fⁱ(x)

                debug && println("0 ∈ fⁱ(x)")

                if 0 ∉ f(interval(prevfloat(X.lo), nextfloat(X.lo)))
                    expansion_pt = X.lo

                elseif 0 ∉ f(interval(prevfloat(X.hi), nextfloat(X.hi)))
                    expansion_pt = X.hi

                else
                    x1 = mid(interval(X.lo, mid(X)))
                    x2 = mid(interval(mid(X), X.hi))
                    if 0 ∉ f(interval(prevfloat(x1), nextfloat(x1))) || 0 ∉ f(interval(prevfloat(x2), nextfloat(x2)))
                        push!(L, interval(X.lo, m))
                        push!(L, interval(m, X.hi))
                        continue

                    else
                        push!(R, Root(X, :unknown))

                        debugroot && @show "Multiple root found", X

                        continue
                    end
                end

            else
                # 0 ∉ fⁱ(x)

                debug && println("0 ∉ fⁱ(x)")

                if (diam(X)/mag(X)) < reltol && diam(f(X)) < abstol
                    push!(R, Root(X, :unknown))

                    debugroot && @show "Tolerance root found", X

                    continue
                end
            end
            # Step 8

            if isinf(expansion_pt)
                expansion_pt = mid(X)
            end


            a = f(interval(expansion_pt))
            b = f′(X)

            if 0 < b.hi && 0 > b.lo && 0 ∉ a
                if a.hi < 0
                    push!(L, X ∩ (expansion_pt - interval(-Inf, a.hi / b.hi)))
                    push!(L, X ∩ (expansion_pt - interval(a.hi / b.lo, Inf)))

                elseif a.lo > 0
                    push!(L, X ∩ (expansion_pt - interval(-Inf, a.lo / b.lo)))
                    push!(L, X ∩ (expansion_pt - interval(a.lo / b.hi, Inf)))

                end

                continue

            else
                N = expansion_pt - (f(interval(expansion_pt))/f′(X))

                debug && (print("Newton step2: "); @show (X, X ∩ N))

                X = X ∩ N
                m = mid(X)

                if isempty(X)
                    continue
                end
            end

            push!(L, X)
        end
    end

    return R
end


"""`newton1d` performs the interval Newton method on the given function `f` and initial interval `x`.
Optional keyword arguments give the tolerances `reltol` and `abstol`.
`reltol` is the tolerance on the relative error whereas `abstol` is the tolerance on |f(X)|,
and a `debug` boolean argument that prints out diagnostic information."""

newton1d{T}(f::Function, x::Interval{T};  args...) =
    newton1d(f, x->D(f,x), x; args...)
