"""
Preconditions the matrix A and b with the inverse of mid(A)
"""
function preconditioner(A::AbstractMatrix, b::AbstractArray)

    Aᶜ = mid.(A)
    B = inv(Aᶜ)

    return B*A, B*b

end

function gauss_seidel_interval(A::AbstractMatrix, b::AbstractArray; precondition=true, maxiter=100)

    n = size(A, 1)
    x = similar(b)
    x .= -1e16..1e16
    gauss_seidel_interval!(x, A, b, precondition=precondition, maxiter=maxiter)
    return x
end
"""
Iteratively solves the system of interval linear
equations and returns the solution set. Uses the
Gauss-Seidel method (Hansen-Sengupta version) to solve the system.
Keyword `precondition` to turn preconditioning off.
Eldon Hansen and G. William Walster : Global Optimization Using Interval Analysis - Chapter 5 - Page 115
"""
function gauss_seidel_interval!(x::AbstractArray, A::AbstractMatrix, b::AbstractArray; precondition=true, maxiter=100)

    precondition && ((A, b) = preconditioner(A, b))

    n = size(A, 1)

    @inbounds for iter in 1:maxiter
        x¹ = copy(x)
        for i in 1:n
            Y = b[i]
            for j in 1:n
                (i == j) || (Y -= A[i, j] * x[j])
            end
            Z = extended_div(Y, A[i, i])
            x[i] = hull((x[i] ∩ Z[1]), x[i] ∩ Z[2])
        end
        if all(x .== x¹)
            break
        end
    end
    x
end

function gauss_seidel_contractor(A::AbstractMatrix, b::AbstractArray; precondition=true, maxiter=100)

    n = size(A, 1)
    x = similar(b)
    x .= -1e16..1e16
    x = gauss_seidel_contractor!(x, A, b, precondition=precondition, maxiter=maxiter)
    return x
end

function gauss_seidel_contractor!(x::AbstractArray, A::AbstractMatrix, b::AbstractArray; precondition=true, maxiter=100)

    precondition && ((A, b) = preconditioner(A, b))

    n = size(A, 1)

    diagA = Diagonal(A)
    extdiagA = copy(A)
    for i in 1:n
        if (typeof(b) <: SArray)
            extdiagA = setindex(extdiagA, Interval(0), i, i)
        else
            extdiagA[i, i] = Interval(0)
        end
    end
    inv_diagA = inv(diagA)

    for iter in 1:maxiter
        x¹ = copy(x)
        x = x .∩ (inv_diagA * (b - extdiagA * x))
        if all(x .== x¹)
            break
        end
    end
    x
end

function gauss_elimination_interval(A::AbstractMatrix, b::AbstractArray; precondition=true)

    n = size(A, 1)
    x = similar(b)
    x .= -1e16..1e16
    x = gauss_elimination_interval!(x, A, b, precondition=precondition)

    return x
end
"""
Solves the system of linear equations using Gaussian Elimination,
with (or without) preconditioning. (kwarg - `precondition`)
Luc Jaulin, Michel Kieffer, Olivier Didrit and Eric Walter - Applied Interval Analysis - Page 72
"""
function gauss_elimination_interval!(x::AbstractArray, a::AbstractMatrix, b::AbstractArray; precondition=true)

    if precondition
        ((a, b) = preconditioner(a, b))
    else
        a = copy(a)
        b = copy(b)
    end
    n = size(a, 1)

    p = zeros(shape(x))

    for i in 1:(n-1)
        if 0 ∈ a[i, i] # diagonal matrix is not invertible
            p .= entireinterval(b[1])
            return p .∩ x
        end

        for j in (i+1):n
            α = a[j, i] / a[i, i]
            b[j] -= α * b[i]

            for k in (i+1):n
                a[j, k] -= α * a[i, k]
            end
        end
    end

    for i in n:-1:1
        sum = 0
        for j in (i+1):n
            sum += a[i, j] * p[j]
        end
        p[i] = (b[i] - sum) / a[i, i]
    end

    p .∩ x
end

function gauss_elimination_interval1(A::AbstractMatrix, b::AbstractArray; precondition=true)

    n = size(A, 1)
    x = fill(-1e16..1e16, n)

    x = gauss_elimination_interval1!(x, A, b, precondition=precondition)

    return x
end
"""
Using `Base.\``
"""
function gauss_elimination_interval1!(x::AbstractArray, a::AbstractMatrix, b::AbstractArray; precondition=true)

    precondition && ((a, b) = preconditioner(a, b))

    a \ b
end

\(A::StaticMatrix{Interval{T}}, b::StaticArray{Interval{T}}; kwargs...) where T = gauss_elimination_interval(A, b, kwargs...)
\(A::Matrix{Interval{T}}, b::Array{Interval{T}}; kwargs...) where T = gauss_elimination_interval(A, b, kwargs...)
