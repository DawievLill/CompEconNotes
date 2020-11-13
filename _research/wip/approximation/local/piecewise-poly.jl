
# Piecewise polynomial interpolation

using Interpolations
using LinearAlgebra
using Plots
using Polynomials

n = 12
t = LinRange(-1,1,n+1)
y = @. t^2 + t + 0.5*sin(20*t)

scatter(t,y,label="data",leg=:top)

p = LinearInterpolation(t,y) # One option is to have linear interpolation 
plot!(x->p(x),-1,1,label="piecewise linear")

p = CubicSplineInterpolation(t,y) # Alternatively one could have a interpolant that is piecewise cubic
plot!(x->p(x),-1,1,label="piecewise cubic")