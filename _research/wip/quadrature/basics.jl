
# Basic integration techniques in Julia

using Plots
using QuadGK

exact = exp(1) - 1 # provides the exact solution to the problem. 

Q, errest = quadgk(x -> exp(x), 0, 1)
@show Q # The result from this package is exactly the same as the exact solution

# However, numerical integration can easily be applied to function that do not provide analytical solutions. 

Q,errest = quadgk(x -> exp(sin(x)), 0, 1)
@show Q

plot([exp,x->exp(sin(x))],0,1,fill=0,leg=:none, 
    ylabel=["exp(x)" "exp(sin(x))"],ylim=[0,2.7],layout=(2,1))
