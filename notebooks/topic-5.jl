### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ eeca0160-d37e-11eb-05f4-bd7469319dad
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
			Pkg.PackageSpec(name="Images"), 
			Pkg.PackageSpec(name="ImageMagick"), 
			Pkg.PackageSpec(name="PlutoUI"), 
			Pkg.PackageSpec(name="HypertextLiteral"),
			Pkg.PackageSpec(name="Random"), 
			Pkg.PackageSpec(name="BenchmarkTools"), 
			Pkg.PackageSpec(name="Plots"), 
			Pkg.PackageSpec(name="LaTeXStrings")
			])
	using Random
	using BenchmarkTools
	using Images
	using PlutoUI
	using HypertextLiteral
	using LinearAlgebra
	using Plots
	using LaTeXStrings
end

# ╔═╡ 99e0b1b9-d032-41d6-8a1d-742050d5ddeb
html"""
<style>
  main {
    max-width: 900px;
  }
</style>
"""

# ╔═╡ 0ac4cb74-c027-4677-bd45-672f19035ce4
md" # Function approximation "

# ╔═╡ 25b6e79e-37dc-4a67-8b56-fbaf9da41713
md" In many of our economic applications we are often confronted with functions that don't have an analytic representation or are difficult to compute. This will be true for the value functions in our dynamic programming discussion for the next session. In order to numerically represent problematic functions on our computer we use approximations, which utilises a simpler function to approximate the more difficult one. Approximations will use data on the function in question. Normally we know some of the function values that correspond to some input points. 

The **goal** is to take our approximation to the data and determine what the function value is on some point that is outside of the finite set of points that we originally evaluated the function on. In economics we do this all the time. We make predictions based on observed data. The big difference with approximation is that we can choose the number of points at which we want to evaluate the function. We produce the data in function approximation. "

# ╔═╡ a4f24d72-fb51-47cc-827d-c9f3cb367998
md" Following the work of Judd, we will focus on three general approaches to approximation. 

1. Local approximations (perturbation methods)
2. Interpolation or collocation
3. Regression

In the case of local approximations we can approximate the function and it's derivative at a single point using Taylor series expansions. With interpolation we require that an approximation must pass through certain points of the function. The `residual` at each grid point needs to be zero with this class of methods. Finally, in the case of regression we would like to minimise some notion of distance without requiring that the points pass through the function in question. "

# ╔═╡ 9ffc49e0-056b-4029-87ce-2fee216277a5
md"  ## Basis functions "

# ╔═╡ 31bc8646-bf07-4a22-9030-196af852251c
md" Before we start with the idea of a basis function, we need to understand what a basis is in the context of linear algebra. For a good description of linear algebra concepts with code, you can visit the [QuantEcon](https://python.quantecon.org/linear_algebra.html) page. We start with the idea of span and the move to linear independence and finally talk about the definition of a basis.

If we are provided a set of vectors $A:=\left\{a_{1}, \ldots, a_{k}\right\}$ in $\mathbb{R}^{n}$ we want to know which new vectors can be created by performing **linear** operations. New vectors that are created in this way are referred to as *linear combinations* of $A$. We have that $y \in \mathbb{R}^{n}$ is a linear combination of $A$ if 

$y=\beta_{1} a_{1}+\cdots+\beta_{k} a_{k} \quad \text{for some scalars} \quad \beta_{1}, \ldots, \beta_{k}$

These $\beta$ values are the coefficients of the linear combination. The set of linear combinations of $A$ is referred to as the **span** of $A$

The following figure shows the span of $A=\left\{a_{1}, a_{2}\right\}$ in $\mathbb{R}^{3}$. It provides a 2 dimensional plane passing through these two points and the origin. 
"

# ╔═╡ b6f096f4-9857-4d66-a62d-8dd818c4d9a2
begin
	
	plotly()
	# fixed linear function, to generate a plane
	f(x, y) = 0.2x + 0.1y
	
	# lines to vectors
	x_vec = [0 0; 3 3]
	y_vec = [0 0; 4 -4]
	z_vec = [0 0; f(3, 4) f(3, -4)]
	
	# draw the plane
	n = 20
	grid = range(-5, 5, length = n)
	z2 = [ f(grid[row], grid[col]) for row in 1:n, col in 1:n ]
	plot(grid, grid, z2, fill = :ice, st = :surface, alpha = 0.75)
	plot!(x_vec, y_vec, z_vec, colour = [:black :black], linewidth = 5, labels = "", colorbar = false)
	
end

# ╔═╡ Cell order:
# ╟─eeca0160-d37e-11eb-05f4-bd7469319dad
# ╟─99e0b1b9-d032-41d6-8a1d-742050d5ddeb
# ╟─0ac4cb74-c027-4677-bd45-672f19035ce4
# ╟─25b6e79e-37dc-4a67-8b56-fbaf9da41713
# ╟─a4f24d72-fb51-47cc-827d-c9f3cb367998
# ╟─9ffc49e0-056b-4029-87ce-2fee216277a5
# ╟─31bc8646-bf07-4a22-9030-196af852251c
# ╠═b6f096f4-9857-4d66-a62d-8dd818c4d9a2
