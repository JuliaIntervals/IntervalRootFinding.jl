"""
Complex numbers as 2-vectors, enough for polynomials.
"""
struct Compl{T}
    re::T
    im::T
end

Base.show(io::IO, c::Compl) = println(io, c.re, " + ", c.im, "im")

Base.:+{T}(a::Compl{T}, b::Compl{T}) = Compl{T}(a.re+b.re, a.im+b.im)
Base.:-{T}(a::Compl{T}, b::Compl{T}) = Compl{T}(a.re-b.re, a.im-b.im)

Base.:*{T}(a::Compl{T}, b::Compl{T}) = Compl{T}(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re)
Base.:*{T}(α::Number, z::Compl{T}) = Compl{T}(α*z.re, α*z.im)

Base.one{T}(z::Compl{T}) = Compl{T}(one(T), zero(T))

Base.copy{T}(z::Compl{T}) = Compl{T}(z.re, z.im)

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
