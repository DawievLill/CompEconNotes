
# Poor conditioning in polynomial interpolation

using LinearAlgebra
using Plots
using Polynomials

n = 5
t = LinRange(-1, 1, n + 1)
y = @. t^2 + t + 0.05 * sin(20 * t) # using the broadcast operator across all the different operations. Don't alway like doing this, but here it seems to cut down on typing... Biggest issue is with sin.(20 .* t). 

scatter(t, y, label="data", leg=:top)

## We can use polyfit to get a polynomial interpolant

p = fit(t, y, n)
plot!(p, -1, 1, label="interpolant") # This polynomial fits the points really well...

# However, if we fit a different set of points then the story changes a bit

n = 18
t = LinRange(-1, 1, n + 1)
y = @. t^2 + t + 0.05 * sin(20 * t)

scatter(t, y, m=:o, l=nothing, label="data", leg=:top)

p = fit(t, y, n)
x = LinRange(-1, 1, 1000)  # use a lot of points
plot!(x, p.(x), label="interpolant") # polynomial interpolant is a terrible fit in this case!