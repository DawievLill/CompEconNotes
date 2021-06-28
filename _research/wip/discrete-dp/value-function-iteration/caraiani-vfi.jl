using Base: Float64
using Parameters

## Implementation of the value function algorithm. Follows the work of Collard. Optimal growth model with a CRRA utility function and law of motion provided for capital. 

# This is an easy implementation, perhaps one of the simplest out of all that I have read. There are other methods, this one is just quite easy to understand

## Start by specifying the parameters

@with_kw struct params
    σ::Float64 = 1.5
    δ::Float64 = 0.1
    β::Float64 = 0.95
    α::Float64 = 0.30
    nbk::Float64 = 1000         # Number of grid points
    dev::Float64 = 0.9          # Maximal deviation from the steady state


p = params()

function cara_vfi()
    @unpack σ, δ, β, α, dev, nbk = p

    kstar = ((1 - β * (1 - δ))/(α * β))^(1/(α - 1))     # Steady state value for capital
    csy = (1 - α * β * δ) / (1 - β * (1 - δ))           # 
    kmin = (1 - dev) * kstar  # 
    kmax = (1 + dev) * kstar

    devk = (kmax - kmin) / (nbk - 1)

    k = collect(range(kmin, kmax, step = nbk))      # Alternatively one could write k = collect(kmin:nbk:kmax)
    

