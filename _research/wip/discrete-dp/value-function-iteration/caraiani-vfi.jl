
using Parameters

## Implementation of the value function algorithm. Follows the work of Collard. Optimal growth model with a CRRA utility function and law of motion provided for capital. 

# This is an easy implementation, perhaps one of the simplest out of all that I have read. There are other methods, this one is just quite easy to understand

## Start by specifying the parameters

@with_kw struct params
    
    # parameter specification
    σ       = 1.5            # utility parameter
    δ       = 0.1            # depreciation rate
    β       = 0.95           # discount factor
    α       = 0.30           # capital elasticity of output

    # steady state
    ks      = ((1 - β * (1 - δ)) / (α * β)) ^ (1 / (α - 1))

    # grid construction
    nbk     = 1000           # Number of grid points
    dev     = 0.9            # Maximal deviation from the steady state
    kmin    = (1 - dev) * ks # Lower bound on the grid 
    kmax    = (1 + dev) * ks # Upper bound on the grid

p = params()

function cara_vfi()
    @unpack σ, δ, β, α, dev, nbk, kmin, kmax = p

    csy = (1 - α * β * δ) / (1 - β * (1 - δ))            
    devk = (kmax - kmin) / (nbk - 1)

    k = collect(range(kmin, kmax, step = nbk))      # Alternatively one could write k = collect(kmin:nbk:kmax)

end
    

