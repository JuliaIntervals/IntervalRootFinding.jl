function realify(f)
    function g(x)
        z2 = f(Complex(x...))
        return [z2.re, z2.im]
    end

    return g
end

"""
Takes the derivative of a complex function and returns the real jacobian
that implements it.
"""
function realify_derivative(fp)
    function g_jac(x)
        fpz = fp(Complex(x...))
        return [fpz.re -fpz.im ; fpz.im fpz.re]
    end
    return g_jac
end
