
using LinearAlgebra
using Parameters

# Value function iteration
# This is basically just a port of the MATLAB code. Trying to understand the basic ideas. 
# Setup of the model parameters

@with_kw struct params
    # basic parameters
    α = 1/3             # Capital share in production function
    β = 0.95            # Time discount factor
    δ = 0.05            # Depreciation rate
    σ = 1.              # CRRA (=1/IES)
    ρ = (1/β) - 1       # This is the implied rate of time preference 
    
    # grid construction
    max_iter = 500      # Maximum number of iterations
    tol = 1e-7          # Numbers smaller than this are treated as zero
    penalty = 10^16     # Penalises constraint violations (will have to think about this)
    n = 1001            # Number of nodes for `k` grid
    
    # steady state
    kstar = (α/(ρ .+ δ))^(1/(1 .- α)) # Steady state Value
    kbar  = (1/δ)^(1/(1 .- α))
end

p = params()

function ed_vfi()
    @unpack α, β, δ, σ, ρ, max_iter, tol, penalty, n, kstar, kbar = p      

    # Components for the grid of capital stock
    kmin = tol # Lower bound for `k`
    kmax = kbar # Upper bound for `k`

    k = LinRange(kmin, kmax, n) |> collect  # These area linearly spaced, but we need to choose the grid better / more carefully. 
    
    c = zeros(n, n) # This is a nxn matrix that contains only zeros. 
    
    # The for loop works at this point. Just need to figure out what it is doing. 
    for j in 1:n
        c[:, j] = k.^α .+ (1 .- δ) .* k .- k[j]
    end 

    # This for loop leads to some infeasible choices, so we need to enforce feasibility. We will do this by penalising violations of the feasibility constraint. This section needs some work. 

    violations = [c .<= 0]
    c = c[c .> 0] .+ eps()

    if σ == 1.
        u = log.(c) .- (penalty .* violations[1]);
    else
        u = (1 ./(1 .- σ))*(c.^(1 .- σ) .- 1) .- (penalty .* violations);
    end

    # Solve the Bellman equation by value function iteration
    v = zeros(n, 1)

    for i = 1:max_iter
        RHS = u .+ β .* kron(ones(n, 1), v') 



end