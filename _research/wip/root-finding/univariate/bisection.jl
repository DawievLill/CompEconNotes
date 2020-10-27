
#=
Bisection algorithm

@author: Dawie van Lill <dvanlill@sun.ac.za>

@date: 2020-10-27

References
----------
Wheeler - Algorithms for Optimisation p49
Miranda and Fackler -- Applied Computational Economics and Finance

=#

using Parameters

@with_kw struct params
    tol::Float64 = 1e-8
end

f = x -> x^3 - 2;
p = params()

# Bisection algorithm as presented in Miranda and Fackler
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

## Example 
fackler_bisection(f, 1, 2) # This should give a value of about 1.25992105 

# This algorithm gives us a bracket that contains the root
function bisection_wheeler(f, a, b)
    @unpack tol = p

    # Interesting way in which an if statement can be phrased in one line
    if a > b; a, b = b, a; end # This ensures that a < b

    ya, yb = f(a), f(b)

    if ya == 0; b = a; end       
    if yb == 0; a = b; end

    while b - a > tol
        x = (a + b) / 2
        y = f(x)

        if y == 0
            a, b = x, x
        elseif sign(y) == sign(ya)
            a = x 
        else 
            b = x
        end

    end
    return (a, b)   
end



