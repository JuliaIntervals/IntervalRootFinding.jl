"""
Complex numbers as 2-vectors, enough for polynomials.
"""
struct Compl{T} <: Number
    re::T
    im::T
end

Base.promote_rule(::Type{Compl{T}}, ::Type{S}) where {T, S<:Real} = Compl{T}
Base.convert(::Type{Compl{T}}, x::Real) where {T} = Compl{T}(x, zero(T))
# Base.convert(::Type{Compl{T}}, x::Real) where {S,T<:Interval{S}} = Compl{T}(interval(S, x), zero(T))


Base.show(io::IO, c::Compl) = println(io, c.re, " + ", c.im, "im")

Base.:+(a::Compl{T}, b::Compl{T}) where {T} = Compl{T}(a.re+b.re, a.im+b.im)
Base.:-(a::Compl{T}, b::Compl{T}) where {T} = Compl{T}(a.re-b.re, a.im-b.im)

Base.:*(a::Compl{T}, b::Compl{T}) where {T} = Compl{T}(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re)


Base.:*(α::Number, z::Compl{T}) where {T} = Compl{T}(α*z.re, α*z.im)

Base.one(z::Compl{T}) where {T} = Compl{T}(one(T), zero(T))

Base.copy(z::Compl{T}) where {T} = Compl{T}(z.re, z.im)

"""
Takes a complex (polynomial) function f and returns a function g:R^2 -> R^2
that implements it.
"""
function realify(f)

    function g(x)
        z = Compl(x[1], x[2])
        z2 = f(z)
        SVector(z2.re, z2.im)
    end

    return g

end

"""
Takes the derivative of a complex function and returns the real jacobian
that implements it.
"""
function realify_derivative(fp)
    function g_jac(x)
        z = Compl(x[1], x[2])
        fpz = fp(z)
        SMatrix{2, 2}(fpz.re, fpz.im, -fpz.im, fpz.re)
    end
    return g_jac
end
