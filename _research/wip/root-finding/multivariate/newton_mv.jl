#=
Newton's method (multivariate case)

@author: Dawie van Lill <dvanlill@sun.ac.za>

@date: 2020-10-27

References
----------
Fundamentals of numerical computation

=#

using LinearAlgebra
using Parameters

@with_kw struct params
    tol::Float64     = 1000 * eps()
    xtol::Float64    = 1000 * eps()
    maxiter::Int64   = 40
end

p = params()

"""
newtonsys(f,jac,x1)

Use Newton's method to find a root of a system of equations,
starting from `x1`. The functions `f` and `jac should return the
residual vector and the Jacobian matrix, respectively. Returns
history of root estimates as a vector of vectors.
"""
function newtonsys(f, jac, x1)
    @unpack tol, xtol, maxiter = p

    x = [float(x1)]     # Convert to float since the eventual output will most likely be floating point
    y, J = f(x1), jac(x1)
    dx = Inf
    k = 1

    # The norm is used instead of the absolute value
    while (norm(dx) > xtol) && (norm(y) > tol) && (k < maxiter)

        dx = -(J\y)     # Newton step uses backslash
        push!(x, x[k] + dx)
        k += 1
        y, J = f(x[k]), jac(x[k])
    end
    return x
end

# Looks incredibly similar to the univariate case. This is a nice feature of Julia, the fact that it is vector-friendly. 

# Example

function nlfun(x)
    f = zeros(3)  
    f[1] = exp(x[2]-x[1]) - 2;
    f[2] = x[1]*x[2] + x[3];
    f[3] = x[2]*x[3] + x[1]^2 - x[2];
    return f
end
   
# Jacobian is manually calculated here. There should be a way to do this automatically. 
function nljac(x)
    J = zeros(3,3)
    J[1,:] = [-exp(x[2]-x[1]),exp(x[2]-x[1]), 0]
    J[2,:] = [x[2], x[1], 1]
    J[3,:] = [2*x[1], x[3]-1, x[2]]
    return J
end

x1 = [0,0,0]
x = newtonsys(nlfun,nljac,x1)