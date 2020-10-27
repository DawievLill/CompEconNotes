
#=
Bisection algorithm

@author: Dawie van Lill <dvanlill@sun.ac.za>

@date: 2020-10-27

References
----------
Wheeler - Algorithms for Optimisation p49
Miranda and Fackler -- Applied Computational Economics and Finance
Driscoll -- Fundamentals of Numerical Computation

=#

using Parameters

@with_kw struct params
    tol::Float64 = 1e-8
end

f = x -> x^3 - 2;
p = params()

function fackler_bisection(f, a, b)
    @unpack tol = p

    s = sign(f(a))
    x = (a + b) / 2
    d = (b - a) / 2

    while d > tol
        d = d/2

        if s == sign(f(x))
            x = x + d
        else
            x = x - d
        end
    end
    return x
end





