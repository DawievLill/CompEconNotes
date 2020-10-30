#=
Newtons method for scalar root finding

@author: Dawie van Lill <dvanlill@sun.ac.za>

@date: 2020-10-27

References
----------
Fundamentals of numerical computation

=#

using ForwardDiff 
using Parameters
using Roots

@with_kw struct params
    funtol::Float64 = 100 * eps()
    xtol::Float64   = 100 * eps()
    maxiter::Int64  = 40
end

p = params()

"""
newton(f, dfdx, x1)

Use Newton's method to find a root of `f` starting from `x1`, where `dfdx` is the derivative of `f`. This returns a vector of root estimates.
"""
function newton(f, dfdx, x1)
    @unpack funtol, xtol, maxiter = p

    x  = [x1]    # This puts the `x1` value into an array
    y  = f(x1)   # This is for a specifically defined input function
    dx = Inf     # For the initial pass -- think about This
    k = 1

    while (abs(dx) > xtol) && (abs(y) > funtol) && (k < maxiter)

        dydx = dfdx(x[k])
        dx   = -y/dydx      # This is referred to as the Newton step
        push!(x, x[k] + dx) # Amend this array by adding the adjusted Newton step value 

        k += 1              # Increment k by one (next iteration)
        y = f(x[k])         # Function evaluation at the next value of the iteration

    end

    if k == maxiter
        @warn "Maximum number of iterations reached"
    end
    return x
end

# Example 

f(x) = exp(x) - x - 2;        # Our function that we want to evaluate
dfdx(x) = exp(x) - 1;         # We can also determine this with AD of symbolic differentiation. 
x = newton(f, dfdx, 1.0);        
# println(x[end])               # Provide us with the value for the root


# There is an implementation of Newton's method using the Roots.jl package and also ForwardDiff to find the derivative

D(f) = x -> ForwardDiff.derivative(f, float(x))
@time find_zero((f, D(f)), 1.0, Roots.Newton())