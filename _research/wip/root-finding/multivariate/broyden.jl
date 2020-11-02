#=
Broyden's method

@author: Dawie van Lill <dvanlill@sun.ac.za>

@date: 2020-10-27

References
----------
Fundamentals of numerical computation
https://nicoguaro.github.io/posts/numerical-05/#disqus_thread


Should be quite close to Newton's method
=#

using LinearAlgebra
using Parameters

@with_kw struct params
    maxiter::Int64 = 50
    tol::Float64   = 1000 * eps
    xtol::Float64  = 1000 * eps 
end

p = params()

"""
broyden(f, jac, x)

Use Broyden's method to find a root of a system of equations,
starting from `x`. The functions `f` and `jac' should return the
residual vector and the Jacobian matrix, respectively. Returns
history of root estimates as a vector of vectors.
"""
function broyden(f, jac, x)
    @unpack maxiter, tol = p

    y, J = f(x), jac(x)
    dx   = Inf
    k    = 1

    # The norm is used instead of the absolute value
    while (norm(dx) > xtol) && (norm(y) > tol) && (k < maxiter)

        dx = -J\y
        push!(x, x[k] + dx)
        k += 1
        y = f(x[k])

        # Need to figure this step out
        df = y .- 
        J = jac(x[k]) .+ (df .- jac(x[k]) .* dx) .* dx' / (dx' .* dx)
    end
    return x
end


function f(x)
    return [x[1] + 2*x[2] - 2, x[1]^2 + 4*x[2]^2 - 4]
end

function jac(x)
    return [1 2;
           2*x[1] 8*x[2]]
end

print(broyden(fun, jac, [1.0, 2.0]))