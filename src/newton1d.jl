doc"""`newton1d` performs the interval Newton method on the given function `f`
with its derivative `f′` and initial interval `x`.
Optional keyword arguments give the tolerances `epsX` and `epsf`.
`epsX` is the tolerance on the relative error whereas `epsf` is the tolerance on |f(X)|,
and a `debug` boolean argument that prints out diagnostic information."""

function newton1d{T}(f::Function, f′::Function, x::Interval{T};
                    epsX=eps(T), epsf=eps(T), debug=false)

    L = Array{Interval, 1}()

    R = Array{Root, 1}()

    push!(L, x)

    while !isempty(L)
        X = pop!(L)
        m = mid(X)

        if 0 ∉ f′(X)
            # O ∉ f'(X)
            while true
                m = mid(X)
                N = m - (f(Interval(m)) / f′(X))
                X = X ∩ N

                if isempty(X)
                    break

                elseif 0 ∈ f(Interval(prevfloat(m), nextfloat(m)))
                    push!(R, Root(X, :unique))
                    break
                end
            end

        else
            # 0 ∈ f'(X)
            expansion_pt = Inf
            # expansion point for the newton step might be m, X.lo or X.hi according to some conditions

            if 0 ∈ f(Interval(mid(X)))
                # 0 ∈ fⁱ(x)
                # Step 7

                if 0 ∉ f(Interval(X.lo))
                    expansion_pt = X.lo

                elseif 0 ∉ f(Interval(X.hi))
                    expansion_pt = X.hi

                else
                    x1 = 0.25 * (3 * X.lo + X.hi)
                    x2 = 0.25 * (X.lo + 3 * X.hi)
                    if 0 ∉ f(Interval(x1)) || 0 ∉ f(Interval(x2))
                        push!(L, Interval(X.lo, m))
                        push!(L, Interval(m, X.hi))
                        break

                    else
                        push!(R, Root(X, :unique))
                        break
                    end
                end

            else
                # 0 ∉ fⁱ(x)

                if (diam(X)/mag(X)) < epsX && diam(f(X)) < epsf
                    push!(R, Root(X, :unknown))
                    break
                end
            end
            # Step 8

            if isinf(expansion_pt)
                expansion_pt = mid(X)
            end

            initial_width = diam(X)

            N = expansion_pt - (f(Interval(expansion_pt))/f′(X))
            X = X ∩ N
            m = mid(X)

            if isempty(X)
                break
            end
            # How can a Newton step in (9.2.2) return two intervals?
            if diam(X) > initial_width/2
                push!(L, Interval(m, X.hi))
                X = Interval(X.lo, m)
            end

            push!(L, X)
        end
    end

    return R
end
