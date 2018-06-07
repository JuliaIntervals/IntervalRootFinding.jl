"""
Preconditions the matrix A and b with the inverse of mid(A)
"""
function preconditioner(A::AbstractMatrix, b::AbstractArray)

    Aᶜ = mid.(A)
    B = pinv(Aᶜ)     # TODO If Aᶜ is singular

    return B*A, B*b

end

function newtonnd(f::Function, deriv::Function, x::IntervalBox{NUM, T}; reltol=eps(T), abstol=eps(T), debug=false, debugroot=false) where {T<:AbstractFloat} where {NUM}    # TODO Incorporate Hull and Box consistencies

    L = IntervalBox{NUM, T}[] # Array to hold the interval boxes still to be processed

    R = Root{IntervalBox{NUM, T}}[] # Array to hold the `root` objects obtained

    push!(L, x) # Initialize
    n = size(X, 1)
    while !isempty(L)   # Until all interval boxes have been processed
        Xᴵ = pop!(L)     # Process next interval box
        if !all(0 .∈ f(Xᴵ))
            continue
        end
        Xᴵ¹ = copy(Xᴵ)
        debug && (print("Current interval popped: ");  println(Xᴵ))

        if (isempty(Xᴵ))
            continue
        end
        if diam(Xᴵ) < reltol
            if max((abs.(IntervalBox(f(Xᴵ))))...) < abstol

                debugroot && @show "Tolerance root found", Xᴵ

                push!(R, Root(Xᴵ, :unknown))
                continue
            else
                continue
            end
        end

        next_iter = false

        while true

            use_B = false
            next_iter = false

            initial_width = diam(Xᴵ)
            debug && (print("Current interval popped: ");  println(Xᴵ))
            if use_B
                # TODO Compute X using B in Step 19
            else
                X = mid(Xᴵ)
            end

            J = SMatrix{n}{n}(deriv(Xᴵ))   # either jacobian or calculated using slopes

            # Xᴵ = IntervalBox((X + linear_hull(J, -f(interval.(X)))) .∩ Xᴵ)
            Xᴵ = IntervalBox((X + (J \ -f(interval.(X)))) .∩ Xᴵ)

            if (isempty(Xᴵ))
                next_iter = true
                break
            end

            if diam(Xᴵ) < reltol
                if max((abs.(IntervalBox(f(Xᴵ))))...) < abstol

                    debugroot && @show "Tolerance root found", Xᴵ
                    next_iter = true
                    push!(R, Root(Xᴵ, :unknown))
                    break
                else
                    next_iter = true
                    break
                end
            end

            if all(0 .∈ f(interval.(X))) && initial_width > 0.9diam(Xᴵ)
                next_iter = true
                push!(R, Root(Xᴵ, :unique))
                break
            end

            if diam(Xᴵ) > initial_width / 8
                break
            end

            #Criterion C


        end

        if next_iter
            continue
        end

        if 0.25 * diam(Xᴵ¹) ≤ max((diam.(Xᴵ¹) .- diam.(Xᴵ))...)
            push!(L, Xᴵ)
        else
            push!(L, bisect(Xᴵ)...)
        end
    end

    R

end

newtonnd(f::Function, x::IntervalBox{NUM, T};  args...)  where {T<:AbstractFloat} where {NUM} = newtonnd(f, x->ForwardDiff.jacobian(f,x), x; args...)
