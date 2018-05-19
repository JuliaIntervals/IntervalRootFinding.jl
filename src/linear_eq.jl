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
    x = fill(-1e16..1e16, n)

    if (typeof(b) <: SArray || typeof(b) <: MArray)
        x = MVector{n}(x)
    end

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
        for i in 1:n
            Y = b[i]
            for j in 1:n
                (i == j) || (Y -= A[i, j] * x[j])
            end
            Z = extended_div(Y, A[i, i])
            x[i] = hull((x[i] ∩ Z[1]), x[i] ∩ Z[2])
        end
    end
    x
end
#
# function preconditioner_static{T, N}(A::MMatrix{N, N, Interval{T}}, b::MVector{N, Interval{T}})
#
#     Aᶜ = mid.(A)
#     B = inv(Aᶜ)
#
#     return B*A, B*b
#
# end
#
#
# function gauss_seidel_interval_static{T, N}(A::MMatrix{N, N, Interval{T}}, b::MVector{N, Interval{T}}; precondition=true, maxiter=100)
#
#     n = size(A, 1)
#     x = @MVector fill(-1e16..1e16, n)
#     gauss_seidel_interval_static!(x, A, b, precondition=precondition, maxiter=maxiter)
#     return x
# end
#
# """
# Iteratively solves the system of interval linear
# equations and returns the solution set. Uses the
# Gauss-Seidel method (Hansen-Sengupta version) to solve the system.
# Keyword `precondition` to turn preconditioning off.
# Eldon Hansen and G. William Walster : Global Optimization Using Interval Analysis - Chapter 5 - Page 115
# """
# function gauss_seidel_interval_static!{T, N}(x::MVector{N, Interval{T}}, A::MMatrix{N, N, Interval{T}}, b::MVector{N, Interval{T}}; precondition=true, maxiter=100)
#
#     precondition && ((A, b) = preconditioner_static(A, b))
#
#     n = size(A, 1)
#
#     for iter in 1:maxiter
#         for i in 1:n
#             Y = b[i]
#             for j in 1:n
#                 (i == j) || (Y -= A[i, j] * x[j])
#             end
#             Z = extended_div(Y, A[i, i])
#             x[i] = hull((x[i] ∩ Z[1]), x[i] ∩ Z[2])
#         end
#     end
#     x
# end
#
# function preconditioner_static1{T, N}(A::SMatrix{N, N, Interval{T}}, b::SVector{N, Interval{T}})
#
#     Aᶜ = mid.(A)
#     B = inv(Aᶜ)
#
#     return B*A, B*b
#
# end
#
#
# function gauss_seidel_interval_static1{T, N}(A::SMatrix{N, N, Interval{T}}, b::SVector{N, Interval{T}}; precondition=true, maxiter=100)
#
#     n = size(A, 1)
#     x = @MVector fill(-1e16..1e16, n)
#     gauss_seidel_interval_static1!(x, A, b, precondition=precondition, maxiter=maxiter)
#     return x
# end
#
# """
# Iteratively solves the system of interval linear
# equations and returns the solution set. Uses the
# Gauss-Seidel method (Hansen-Sengupta version) to solve the system.
# Keyword `precondition` to turn preconditioning off.
# Eldon Hansen and G. William Walster : Global Optimization Using Interval Analysis - Chapter 5 - Page 115
# """
# function gauss_seidel_interval_static1!{T, N}(x::MVector{N, Interval{T}}, A::SMatrix{N, N, Interval{T}}, b::SVector{N, Interval{T}}; precondition=true, maxiter=100)
#
#     precondition && ((A, b) = preconditioner_static1(A, b))
#
#     n = size(A, 1)
#
#     for iter in 1:maxiter
#         for i in 1:n
#             Y = b[i]
#             for j in 1:n
#                 (i == j) || (Y -= A[i, j] * x[j])
#             end
#             Z = extended_div(Y, A[i, i])
#             x[i] = hull((x[i] ∩ Z[1]), x[i] ∩ Z[2])
#         end
#     end
#     x
# end
#
# function preconditioner_static2{T, N}(A::SMatrix{N, N, Interval{T}}, b::SVector{N, Interval{T}})
#
#     Aᶜ = mid.(A)
#     B = inv(Aᶜ)
#
#     return B*A, B*b
#
# end
#
#
# function gauss_seidel_interval_static2{T, N}(A::SMatrix{N, N, Interval{T}}, b::SVector{N, Interval{T}}; precondition=true, maxiter=100)
#
#     n = size(A, 1)
#     x = @SVector fill(-1e16..1e16, n)
#     x = gauss_seidel_interval_static2!(x, A, b, precondition=precondition, maxiter=maxiter)
#     return x
# end
#
# """
# Iteratively solves the system of interval linear
# equations and returns the solution set. Uses the
# Gauss-Seidel method (Hansen-Sengupta version) to solve the system.
# Keyword `precondition` to turn preconditioning off.
# Eldon Hansen and G. William Walster : Global Optimization Using Interval Analysis - Chapter 5 - Page 115
# """
# function gauss_seidel_interval_static2!{T, N}(x::SVector{N, Interval{T}}, A::SMatrix{N, N, Interval{T}}, b::SVector{N, Interval{T}}; precondition=true, maxiter=100)
#
#     precondition && ((A, b) = preconditioner_static2(A, b))
#
#     n = size(A, 1)
#
#     for iter in 1:maxiter
#         for i in 1:n
#             Y = b[i]
#             for j in 1:n
#                 (i == j) || (Y -= A[i, j] * x[j])
#             end
#             Z = extended_div(Y, A[i, i])
#             x = setindex(x, hull((x[i] ∩ Z[1]), x[i] ∩ Z[2]), i)
#         end
#     end
#     x
# end

function gauss_seidel_contractor(A::AbstractMatrix, b::AbstractArray; precondition=true, maxiter=100)

    n = size(A, 1)
    x = fill(-1e16..1e16, n)

    if (typeof(b) <: SArray || typeof(b) <: MArray)
        x = MVector{n}(x)
    end

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
        x = x .∩ (inv_diagA * (b - extdiagA * x))
    end
    x
end