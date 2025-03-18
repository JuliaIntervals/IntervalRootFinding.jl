"""
    preconditioner(A, b)

Preconditions the matrix `A` and `b` with the inverse of `mid(A)`
"""
function preconditioner(A::AbstractMatrix, b::AbstractArray)
    Aᶜ = mid.(A)
    det(Aᶜ) == 0 && return A, b
    B = interval.(inv(Aᶜ))
    return B*A, B*b
end

function gauss_seidel_interval(A::AbstractMatrix, b::AbstractArray; precondition=true, maxiter=100)
    x = similar(b)
    x .= -1e16..1e16
    gauss_seidel_interval!(x, A, b, precondition=precondition, maxiter=maxiter)
    return x
end
"""
    gauss_seidel_interval!(x, A, b; precondition=true, maxiter=100)

Iteratively solves the system of interval linear
equations and returns the solution set. Uses the
Gauss-Seidel method (Hansen-Sengupta version) to solve the system.
Keyword `precondition` to turn preconditioning off.

Reference: [Hansen and Walster (2003), Chapter 5, Page 115](@cite hansen2003)
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

function gauss_elimination_interval(A0::AbstractMatrix, b0::AbstractArray ; precondition=true)
    if precondition
        A0, b0 = preconditioner(A0, b0)
    end

    A = similar(A0)
    A .= A0
    b = similar(b0)
    b .= b0
    n = size(A, 1)

    p = similar(b)
    p .= interval(0)

    if any(a -> in_interval(0, a), diag(A))
        p .= interval(-Inf, Inf)
        return p
    end

    for i in 1:(n-1)
        for j in (i+1):n
            α = A[j, i] / A[i, i]
            b[j] -= α * b[i]

            for k in (i+1):n
                A[j, k] -= α * A[i, k]
            end
        end
    end

    for i in n:-1:1
        temp = zero(b[1])

        for j in (i+1):n
            temp += A[i, j] * p[j]
        end

        p[i] = (b[i] - temp) / A[i, i]
    end

    return p
end

function gauss_elimination_interval(A::AbstractMatrix, B::AbstractMatrix ; kwargs...)
    return mapreduce(hcat, eachcol(B)) do b
        return gauss_elimination_interval(A, b ; kwargs...)
    end
end

"""
    gauss_elimination_interval!(x, A, b; precondition=true)

Solves the system of linear equations using Gaussian Elimination.
Preconditioning is used when the `precondition` keyword argument is `true`.

Reference: [Jaulin et al. (2001), p.72](@cite jaulin2001)
"""
function gauss_elimination_interval!(x::AbstractArray, A::AbstractMatrix, b::AbstractArray; precondition=true)
    if precondition
        A, b = preconditioner(A, b)
    end
    A = copy(A)
    b = copy(b)

    n = size(A, 1)

    p = similar(b)
    p .= 0

    any(a -> in_interval(0, a), diag(A)) && return x

    for i in 1:(n-1)
        for j in (i+1):n
            α = A[j, i] / A[i, i]
            b[j] -= α * b[i]

            for k in (i+1):n
                A[j, k] -= α * A[i, k]
            end
        end
    end

    for i in n:-1:1
        temp = zero(b[1])

        for j in (i+1):n
            temp += A[i, j] * p[j]
        end

        p[i] = (b[i] - temp) / A[i, i]
    end

    return intersect_region(p, x)
end

function gauss_elimination_interval1(A::AbstractMatrix, b::AbstractArray; precondition=true)

    n = size(A, 1)
    x = fill(-1e16..1e16, n)

    x = gauss_elimination_interval1!(x, A, b, precondition=precondition)

    return x
end
"""
    gauss_elimination_interval1!(x, a, b; precondition=true)

Using `Base.\\`
"""
function gauss_elimination_interval1!(x::AbstractArray, a::AbstractMatrix, b::AbstractArray; precondition=true)

    precondition && ((a, b) = preconditioner(a, b))

    a \ b
end

\(A::SMatrix{S, S, Interval{T}}, b::SVector{S, Interval{T}}; kwargs...) where {S, T} = gauss_elimination_interval(A, b, kwargs...)
\(A::Matrix{Interval{T}}, b::Vector{Interval{T}}; kwargs...) where T = gauss_elimination_interval(A, b, kwargs...)
