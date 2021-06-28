
using LinearAlgebra
using Parameters

# Value function iteration
# This is basically just a port of the MATLAB code. Trying to understand the basic ideas. 
# Setup of the model parameters

@with_kw struct params
    α = 1/3             # Capital share in production function
    β = 0.95            # Time discount factor
    δ = 0.05            # Depreciation rate
    σ = 1.              # CRRA (=1/IES)
    max_iter = 500      # Maximum number of iterations
    tol = 1e-7          # Numbers smaller than this are treated as zero
    penalty = 10^16     # Penalises constraint violations
    n = 1001            # Number of nodes for `k` grid
    ρ = (1/β) - 1       # This is the implied rate of time preference (is it good idea to use the β here?)
end

p = params()

function ed_vfi()
    @unpack α, β, δ, σ, ρ, max_iter, tol, penalty, n = p      
    
    kstar = (α/(ρ .+ δ))^(1/(1 .- α)) # Steady state Value
    kbar  = (1/δ)^(1/(1 .- α))       

    # Components for the grid of capital stock
    kmin = tol # Lower bound for `k`
    kmax = kbar # Upper bound for `k`

    k = LinRange(kmin, kmax, n) # These area linearly spaced, but we need to choose the grid better / more carefully. 
    
    c = zeros(n, n) # This is a nxn matrix that contains only zeros. 
    
    # The for loop works at this point. Just need to figure out what it is doing. 
    for j in 1:n
        c[:, j] = k.^α .+ (1 .- δ) .* k .- k[j]
    end 

    # This for loop leads to some infeasible choices, so we need to enforce feasibility. We will do this by penalising violations of the feasibility constraint. 

    

end