

# Value function iteration for the cake eating problem. This is on the grid, and does not use interpolation. This is not fitted value function iteration. We will get to the idea of fitted value function iteration at another stage. 

# 1. Start with initial guess values
# 2. At each iteration compute the Bellman operator for
# 3. Repeat until convergence is reached

# The cake is continuous, and the value function is a function of continuous values so we need to discretise the state variable. 

# Compute value and policy function sequentially, point by point. This

# Might be necesarry to compute the value function between grid points -- we would then need to use something like interpolation and function approximation to fully flesh out the function. 

# In this case we can solve on the grid, which avoids the interpolation process. 

@with_kw struct params
    β = 0.96           # discount factor
    Wbar = 10
    ngrid = 50
    grid = LinRange(eps, Wbarm, ngrid) # This is now the grid for both the state and decision space. 
end

p = params()

# This seems like quite a tricky approach. Perhaps look at other examples first. 

function bellman(V0)
    @unpack grid = p


