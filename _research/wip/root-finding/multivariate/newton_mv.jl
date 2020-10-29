#=
Newton's method (multivariate case)

@author: Dawie van Lill <dvanlill@sun.ac.za>

@date: 2020-10-27

References
----------
Fundamentals of numerical computation

=#

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

    x = [float(x1)]     # Convert to float since the eventual output will most likely be floating point


end

# Looks incredibly similar to the univariate case. This is a nice feature of Julia, the fact that it is vector-friendly. 