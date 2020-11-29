
# Value function iteration

# This is basically just a port of the MATLAB code. Trying to understand the basic ideas. 

using LinearAlgebra
using Parameters

# Setup of the model parameters
@with_kw struct params
    α = 1/3             # Capital share in production function
    β = 0.95            # Time discount factor
    δ = 0.05            # Depreciation rate
    σ = 1               # CRRA (=1/IES)
    max_iter = 500      # Maximum number of iterations
    tol = 1e-7          # Numbers smaller than this are treated as zero
    penalty = 10^16     # Penalises constraint violations
    n = 1001            # Number of nodes for `k` grid
end

p = params()

function ed_vfi()
    @unpack α, β, δ, σ, max_iter, tol, penalty, n = p

    ρ = (1/β) - 1        # This is the implied rate of time preference
    
    kstar = (α/(ρ .+ δ))^(1/(1 .- α)) # Steady state Value
    kbar  = (1/δ)^(1/(1 .- α))       

    kmin = tol           # Lower bound for `k`
    kmax = kbar          # Upper bound for `k`

    k = 

end