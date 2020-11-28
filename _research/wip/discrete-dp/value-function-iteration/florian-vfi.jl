
# Example of discrete value function iteration for a deterministic growth model, similar to optimal growth or cake eating problem

# Start by setting parameters in the parameter block

using LinearAlgebra
using Optim
using Parameters
using Roots

@with_kw struct params
    α::Float64 = 0.65
    β::Float64 = 0.95
    grid_max::Int64 = 2         # This is going to be an upper bound on the state variable
    n::Int64 = 150              # Number of grid points
    N_iter::Int64 = 3000        # Number of iterations
    
end


kgrid::Float64 = 1e-2:(grid_max - 1e-2)/(n-1):grid_max  # Represents the equispaced grid. There are other ways to choose the grid points, which we will see later when we talk about function approximation methods, see the notes of Rudik. 