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
                elseif z in f(Interval(prevfloat(mid_x), nextfloat(mid_x)))
                    push!(R, Root(X, :unique))
                    break
                end
            end
        else
            # 0 ∈ f'(X)
        end

    end
    return R
end
