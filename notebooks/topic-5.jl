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


# ╔═╡ Cell order:
# ╟─eeca0160-d37e-11eb-05f4-bd7469319dad
# ╟─99e0b1b9-d032-41d6-8a1d-742050d5ddeb
# ╟─0ac4cb74-c027-4677-bd45-672f19035ce4
# ╟─25b6e79e-37dc-4a67-8b56-fbaf9da41713
# ╟─a4f24d72-fb51-47cc-827d-c9f3cb367998
# ╟─9ffc49e0-056b-4029-87ce-2fee216277a5
# ╠═31bc8646-bf07-4a22-9030-196af852251c
