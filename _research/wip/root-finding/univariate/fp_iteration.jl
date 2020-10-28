#=
Fixed point (function) iteration (univariate case)

@author: Dawie van Lill <dvanlill@sun.ac.za>

@date: 2020-10-27

References
----------
Fundamentals of numerical computation
https://github.com/AEM7130/spring-2020

=#

using LinearAlgebra

# Example 

f = x -> x^(-0.5)

function rudik_iteration(f, guess)

    tol = 1e-3
    x_old = guess

    x = 1e-10
    difference = 1e10

    @time while abs(difference) > tol

        x = f(x_old)
        difference = x - x_old
        x_old = x
    end

    print("Fixed point of f(x) is at $x")

end

## Fix this code, it doesn't seem to run at this stage. 

function caraiani_iteration(f, guess)

    tol = 1e-3
    maxiter = 1000

    for i = 1:maxiter
        if norm(f(guess) - guess) < tol
            return 
        end
    guess = f(guess)
    end
    guess
end

