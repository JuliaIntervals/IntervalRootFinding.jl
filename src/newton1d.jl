"""`newton1d` performs the interval Newton method on the given function `f`
with its derivative `f′` and initial interval `x`.
Optional keyword arguments give the tolerances `reltol` and `abstol`.
`reltol` is the tolerance on the relative error whereas `abstol` is the tolerance on |f(X)|,
and a `debug` boolean argument that prints out diagnostic information."""

function newton1d{T}(f::Function, f′::Function, x::Interval{T};
                    reltol=10eps(T), abstol=eps(T), debug=false)

    L = Interval{T}[]

    R = Root{Interval{T}}[]

    push!(L, x)

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
                N = m - (f(Interval(m)) / f′(X))

                debug && (print("Newton step: "); @show (X, X ∩ N))

                X = X ∩ N

                if isempty(X)
                    break

                elseif 0 ∈ f(Interval(prevfloat(m), nextfloat(m)))
                    push!(R, Root(X, :unique))

                    debug && @show "Root found", X

                    break
                end
            end

        else
            # 0 ∈ f'(X)

            debug && println("0 ∈ f'(X)")

            expansion_pt = Inf
            # expansion point for the newton step might be m, X.lo or X.hi according to some conditions

            if 0 ∈ f(Interval(mid(X)))
                # 0 ∈ fⁱ(x)

                debug && println("0 ∈ fⁱ(x)")

                if 0 ∉ f(Interval(X.lo))
                    expansion_pt = X.lo

                elseif 0 ∉ f(Interval(X.hi))
                    expansion_pt = X.hi

                else
                    x1 = mid(Interval(X.lo, mid(X)))
                    x2 = mid(Interval(mid(X), X.hi))
                    if 0 ∉ f(Interval(x1)) || 0 ∉ f(Interval(x2))
                        push!(L, Interval(X.lo, m))
                        push!(L, Interval(m, X.hi))
                        continue

                    else
                        push!(R, Root(X, :unknown))

                        debug && @show "Multiple root found", X

                        continue
                    end
                end

            else
                # 0 ∉ fⁱ(x)

                debug && println("0 ∉ fⁱ(x)")

                if (diam(X)/mag(X)) < reltol && diam(f(X)) < abstol
                    push!(R, Root(X, :unknown))

                    debug && @show "Tolerance root found", X

                    continue
                end
            end
            # Step 8

            if isinf(expansion_pt)
                expansion_pt = mid(X)
            end

            initial_width = diam(X)

            a = f(Interval(expansion_pt))
            b = f′(X)

            if 0 < b.hi && 0 > b.lo && 0 ∉ a
                if a.hi < 0
                    push!(L, X ∩ (expansion_pt - Interval(-Inf, a.hi / b.hi)))
                    push!(L, X ∩ (expansion_pt - Interval(a.hi / b.lo, Inf)))

                elseif a.lo > 0
                    push!(L, X ∩ (expansion_pt - Interval(-Inf, a.lo / b.lo)))
                    push!(L, X ∩ (expansion_pt - Interval(a.lo / b.hi, Inf)))

                end

                continue

            else
                N = expansion_pt - (f(Interval(expansion_pt))/f′(X))

                debug && (print("Newton step: "); @show (X, X ∩ N))

                X = X ∩ N
                m = mid(X)

                if isempty(X)
                    continue
                end
            end

            if diam(X) > initial_width/2
                push!(L, Interval(m, X.hi))
                X = Interval(X.lo, m)
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
