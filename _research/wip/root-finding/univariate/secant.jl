#=
Quasi-Newton (secant) method for scalar root finding

@author: Dawie van Lill <dvanlill@sun.ac.za>

@date: 2020-10-27

References
----------
Fundamentals of numerical computation
Collard notes

=#

using Parameters

@with_kw struct params
    tol::Float64     = 100 * eps()
    xtol::Float64    = 100 * eps()
    maxiter::Int64   = 40
end

p = params()

"""
secant(f, x1, x2)

Use the secant method to find a root of `f` starting from `x1` and `x2`. Returns a vector of root estimates. 
"""
function secant(f, x1, x2)
    @unpack tol, xtol, maxiter = p

    x = [x1, x2]
    y1 = f(x1)
    y2 = 100
    dx = Inf
    k = 2

    while (abs(dx) > xtol) && (abs(y2) > tol) && (k < maxiter)
        
        y2 = f(x[k])
        dx = -y2 * (x[k] - x[k-1]) / (y2 - y1)  # This is the secant step
        push!(x, x[k] + dx)

        k += 1
        y1 = y2

    end
    return x
end

# The code for the above secant function is almost identical to the newton code. 

# Example 

f = x -> exp(x) - x - 2;        # Our function that we want to evaluate
#dfdx = x -> exp(x) - 1;          # We can also determine this with AD of symbolic differentiation. 
x = secant(f, 1.0, 2.0);        
print(x[end])                   # Provide us with the value for the root

# Another way to potentially approach this problem is the following. However, I don't think the coding is as generic (or elegant as the previous approach)

# function caraiani_secant(f, x_0, x_1)

#     tol = 1e-3
#     err = 1.0
    
#     while err > 0 
        
#         d = (-f(x_1) * (x_1 - x_0))/(f(x_1) - f(x_0))
#         x = x_1 + d

#         err = abs(x - x_1) - tol * (1 + abs(x))

#         x_0 = x_1
#         x_1 = x

#     end
#     return x, f(x)
# end

