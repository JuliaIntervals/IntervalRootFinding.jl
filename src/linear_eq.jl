"""
Preconditions the matrix A and b with the inverse of mid(A)
"""
function preconditioner{T}(A::Matrix{Interval{T}}, b::Array{Interval{T}})

    Aᶜ = mid.(A)
    B = pinv(Aᶜ)

    return B*A, B*b

end


"""
Iteratively solves the system of interval linear
equations and returns the solution set. Uses the
Gauss-Seidel method (Hansen-Sengupta version) to solve the system.
Keyword `precondition` to turn preconditioning off.
"""
function gauss_seidel_interval!{T}(x::Array{Interval{T}}, A::Matrix{Interval{T}}, b::Array{Interval{T}}; precondition=true, maxiter=100)

    precondition && ((M, r) = preconditioner(A, b))

    n = size(A, 1)

    for iter in 1:maxiter
        for i in 1:n
            Y = r[i]
            for j in 1:n
                (i == j) || (Y -= M[i, j] * x[j])
            end
            Y = extended_div(Y, M[i, i])
            x[i] = hull((x[i] ∩ Y[1]), x[i] ∩ Y[2])
        end
    end
    x
end
