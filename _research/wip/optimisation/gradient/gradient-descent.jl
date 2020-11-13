
# Gradient descent method

# Code for this can be found at https://nicoguaro.github.io/posts/numerical-06/ and at https://shuvomoy.github.io/blog/programming/2020/04/10/Implementing_simple_gradient_descent_Julia.html

# The ideas are nicely explained in Wheeler. I am assuming that there will be millions of different ways of coding this. 

using LinearAlgebra
using Parameters
using ProximalOperators
using Optim

@with_kw struct params
    niter::Float64 = 50
    gtol::Float64 = 1e-8
end

p = params()

# Think about how this algorithm works in particular
function grad_descent(f, ∇f, x)
    @unpack niter, gtol = p

    g_old = ∇f(x)
    γ = 0.1 # This is the step size

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


# One nice way to approach this is using Optim.jl


# Put a pin in this last method below, I don't know much about proximal operators. 

## An alternative method, which seems more comprehensive. Try and follow along with the logic. 

# For this problem we will be creating types in Julia 

struct GD_problem{F <: ProximableFunction, A <: AbstractVecOrMat{<:Real}, R <: Real}

    # This is the problem structure that contains information regarding the specific gradient descent problem at hand. 

    f::F        # Our objective function is a Proximable Function types
    x0::A       # The initial condition is of type Vector or Matrix with real entries
    γ::R        # The stepsize is a real value 

end

# Example for the usage of this struct. If you want to solve least squares problem.

A = randn(6, 5)
b = randn(6)
m, n  = size(A)

x0 = randn(n)
f = LeastSquares(A, b)
γ = 1.0

problem  = GD_problem(f, x0, γ)
