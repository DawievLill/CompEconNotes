### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ f4226cfe-ee06-4c72-9615-fc4aedfd045c
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
			])
	using Random
	using BenchmarkTools
	using Images
	using PlutoUI
	using HypertextLiteral
	using LinearAlgebra
	using Plots
end

# ╔═╡ 1f4407d0-9b73-11eb-0e91-cd0de83535aa
md" # Optimisation "

# ╔═╡ 16478345-300b-460b-8198-faccaf7740e9
md" We will probably need about two to three sessions to cover the basic optimisation ideas that are relevant for economics. Optimisation is an incredibly important topic in economics, and many related disciplines (such as statistics and machine learning), so it is worthwhile to spend some time exploring the general themes.

Julia has a well-developed ecosystem for optimisation routines and can comfortably used to meet all your optimisation needs. We will speak to the packages relevant for optimisation later in this session. First, we will cover some basics on optimisation and then move to a discussion of several algorithms relevant for optimisation. 

Unfortunately we will only have time to cover a very small subset of the techniques that are relevant to our field. I strongly recommend that you explore this topic in more detail. Some resources will be provided toward the end of the session. The best book on the topic that utilises Julia code is Algorithms for Optimisation, which can be found [here](https://mykel.kochenderfer.com/textbooks/). "

# ╔═╡ 24d0123f-acc4-4656-b564-d7677bd10cf8
md" ## Basics "

# ╔═╡ 9cf3916d-9436-443e-9ed8-75ce7af69d85
md" For this next session I borrow heavily from the excellent notes of Florian Oswald, which can be found [here](https://floswald.github.io/NumericalMethods/lecture3/). Let us start with a gentle reminder of optimisation theory, so that we understand what our algorithms are trying to achieve. We will obviously only scratch the surface of optimisation theory, this is a subject in its own right.

All the algorithms that are used in this section will rely on some form of an iterative procedure. By this we mean that the algorithms will try and improve the value of some objective function over a succession of steps. The way in which the next step is generated is the distinguishing feature of these algorithms. "

# ╔═╡ c1e1db15-eddb-4535-8d00-c1180a91fdae
md" ### Optimisation packages in Julia "

# ╔═╡ 3b20c100-ff95-4771-957f-634f43d2ccd6
md" Before we move to a discussion on the different algorithms, it might be worthwhile to have a broad overview of the packages that are available in Julia that are intended for optimisation."

# ╔═╡ 38a82b17-13d2-435b-85c8-24c2c619d922
md" ## Algorithms for optimisation"

# ╔═╡ 0170f06a-79ba-4ba3-bfe0-db1e63bc7b74
md" ### Root finding algorithms"

# ╔═╡ 46095078-d68c-4570-b2b5-b8aad3656186
md" The goal of the root finding problem can be stated as follows: 

Given a function $f: \mathbb{R}^{n} \rightarrow \mathbb{R}^n$, find the $x \in \mathbb{R}^{n}$ such that $f(x) = \mathbf{0}$. 

**Example:** If you think about the first model you did in first-year economics, you were solving for the prices and quantities. Clearing markets in a general equilibrium model entails finding the prices such that excess demand functions are zero. "

# ╔═╡ 2d69fac9-ad45-428e-9b78-971eaef81e88
md" #### Bisection method "

# ╔═╡ 0a5f2ec5-0419-4b55-9475-b551e56cf7d5
md" ### Unconstrained optimisation "

# ╔═╡ fbc93dbd-58ab-4666-b866-a72a4a2ec866
md" #### Newton's method "

# ╔═╡ 2d004257-c9f2-40de-9607-d73d668f4537
md" #### Quasi-Newton "

# ╔═╡ f3f3e1ee-092e-4270-96f8-0137a4c9f03a
md" #### EM algorithm "

# ╔═╡ 7c117f42-fa04-4556-b1ae-659680e6e06b
md" ### Constrained optimisation "

# ╔═╡ Cell order:
# ╟─f4226cfe-ee06-4c72-9615-fc4aedfd045c
# ╟─1f4407d0-9b73-11eb-0e91-cd0de83535aa
# ╟─16478345-300b-460b-8198-faccaf7740e9
# ╟─24d0123f-acc4-4656-b564-d7677bd10cf8
# ╟─9cf3916d-9436-443e-9ed8-75ce7af69d85
# ╟─c1e1db15-eddb-4535-8d00-c1180a91fdae
# ╟─3b20c100-ff95-4771-957f-634f43d2ccd6
# ╟─38a82b17-13d2-435b-85c8-24c2c619d922
# ╠═0170f06a-79ba-4ba3-bfe0-db1e63bc7b74
# ╟─46095078-d68c-4570-b2b5-b8aad3656186
# ╟─2d69fac9-ad45-428e-9b78-971eaef81e88
# ╟─0a5f2ec5-0419-4b55-9475-b551e56cf7d5
# ╟─fbc93dbd-58ab-4666-b866-a72a4a2ec866
# ╟─2d004257-c9f2-40de-9607-d73d668f4537
# ╟─f3f3e1ee-092e-4270-96f8-0137a4c9f03a
# ╟─7c117f42-fa04-4556-b1ae-659680e6e06b
