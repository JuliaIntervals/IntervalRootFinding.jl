function newton1d{T}(f::Function, f_prime::Function, x::Interval{T};
                    epsX=eps(T), epsf=eps(T), debug=false)

    L = Interval[]
    R = Root[]
    z = zero(x.lo)
    push!(L, x)
    while !isempty(L)
        X = pop!(L)
        mid_x = mid(X)
        if !(z in f_prime(X))
            # O ∉ f'(X)
            while true
                mid_x = mid(X)
                N = mid_x - (f(Interval(mid_x))/f_prime(X))
                X = X ∩ N
                if isempty(X)
                    break
                elseif z in f(Interval(prevfloat(mid_x), nextfloat(mid_x))) # What if the root is one of the end points of the interval. That way, this condition will never satisfy even thoughthe interval becomes very small. Explicit check for that?
                    push!(R, Root(X, :unique))
                    break
                end
            end
        else
            # 0 ∈ f'(X)
            pointOfExpansion = Inf
            if z in f(Interval(mid(X)))
                # 0 ∈ fⁱ(x)
                # Step 7

                if !(z in f(Interval(X.lo)))
                    pointOfExpansion = X.lo
                elseif !(z in f(Interval(X.hi)))
                    pointOfExpansion = X.hi
                else
                    x1 = 0.25 * (3 * X.lo + X.hi)
                    x2 = 0.25 * (X.lo + 3 * X.hi)
                    if !(z in f(Interval(x1))) || !(z in f(Interval(x2)))
                        push!(L, Interval(X.lo, mid_x))
                        push!(L, Interval(mid_x, X.hi))
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
            if isinf(pointOfExpansion)
                pointOfExpansion = mid(X)
            end
            initial_width = diam(X)
            N = pointOfExpansion - (f(Interval(pointOfExpansion))/f_prime(X))
            X = X ∩ N
            mid_x = mid(X)
            if isempty(X)
                break
            end
            # How can a Newton step in (9.2.2) return two intervals?
            if diam(X) > initial_width/2
                push!(L, Interval(mid_x, X.hi))
                X = Interval(X.lo, mid_x)
            end
            push!(L, X)
        end
    end
    return R
end
