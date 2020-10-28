#=
Bisection algorithm

@author: Dawie van Lill <dvanlill@sun.ac.za>

@date: 2020-10-27

References
----------
Wheeler - Algorithms for Optimisation p49
Miranda and Fackler -- Applied Computational Economics and Finance p33

=#

using BenchmarkTools
using Parameters
using Roots

@with_kw struct params
    tol::Float64 = 1e-18
end

f = x -> x^3 - 2;
p = params()

# Bisection algorithm as presented in Miranda and Fackler
function fackler_bisection(f, a, b)
    @unpack tol = p

    s = sign(f(a))
    x = (a + b) / 2         # Guess for the value of x -- mean of the interval
    d = (b - a) / 2         # Initial difference between bounds

    while d > tol
        d = d/2

        if s == sign(f(x))
            x = x + d
        else
            x = x - d
        end
    end
    print("The root of f(x) is at $x")
end

## Example 
fackler_bisection(f, 1, 2) # This should give a value of about 1.25992105 

# This algorithm gives us a bracket that contains the root
function wheeler_bisection(f, a, b)
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

## Example 
wheeler_bisection(f, 1, 2) # This should give a value of about (1.2599210441112518, 1.2599210515618324)

# There is a package in Julia called roots, that uses something like a bisection algorithm. In fact, you can go look at the code https://github.com/JuliaMath/Roots.jl/blob/master/src/bracketing.jl

# It is worthwhile looking at the code to see how software engineers will code something like this. The general rule is that if there is a package out there to that does things efficiently, do not try and code it up yourself. Use existing software. 

find_zero(f, (1, 2), Bisection()) 

@btime wheeler_bisection(f, 1, 2)
@btime fackler_bisection(f, 1, 2)
@btime find_zero(f, (1, 2), Bisection())

# Comparison of times actually shows that the find_zero method is slowest... It could be that the tolerance is quite high (greater precision).