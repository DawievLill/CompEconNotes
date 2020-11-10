
# Gradient descent method

# Code for this can be found at https://nicoguaro.github.io/posts/numerical-06/ and at https://shuvomoy.github.io/blog/programming/2020/04/10/Implementing_simple_gradient_descent_Julia.html

# The ideas are nicely explained in Wheeler. I am assuming that there will be millions of different ways of coding this. 

using LinearAlgebra
using Parameters
using ProximalOperators

@with_kw struct params
    niter::Float64 = 50
    gtol::Float64 = 1e-8
end

p = params()

# Think about how this algorithm works in particular
function grad_descent(f, ∇f, x)
    @unpack niter, gtol = p

    g_old = ∇f(x)
    γ = 0.1

    for cont = 1:niter

        dx = -γ * g_old
        x = x + dx
        g = ∇f(x)

        dg = g - g_old
        g_old = g

        γ = dx' * dg / (dg' * dg)

        if norm(g) < gtol
            break
        end
    end
    return x, f(x)
end


function rosen(x)
    return (1 - x[1])^2 + 100*(x[2] - x[1]^2)^2
end


function rosen_grad(x)
    return [-2*(1 - x[1]) - 400*x[1]*(x[2] - x[1]^2);
            200*(x[2] - x[1]^2)]
end


println(grad_descent(rosen, rosen_grad, [2.0, 1.0]))


