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
md" We will probably need about two to three sessions to cover the basic optimisation ideas that are relevant for economics. Optimisation is an incredibly important topic in economics, and many related disciplines, so it is worthwhile to spend some time exploring the general themes.

Julia has a well-developed ecosystem for optimisation routines and can comfortably used to meet all your optimisation needs. We will speak to the packages relevant for optimisation later in this session. Our first topic is going to be root finding, then we will move to unconstrained optimisation and end with constrained optimisation."

# ╔═╡ 0170f06a-79ba-4ba3-bfe0-db1e63bc7b74
md" ### Root finding"

# ╔═╡ 46095078-d68c-4570-b2b5-b8aad3656186


# ╔═╡ Cell order:
# ╟─f4226cfe-ee06-4c72-9615-fc4aedfd045c
# ╟─1f4407d0-9b73-11eb-0e91-cd0de83535aa
# ╟─16478345-300b-460b-8198-faccaf7740e9
# ╟─0170f06a-79ba-4ba3-bfe0-db1e63bc7b74
# ╠═46095078-d68c-4570-b2b5-b8aad3656186
