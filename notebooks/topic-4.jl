### A Pluto.jl notebook ###
# v0.14.1

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

# ╔═╡ Cell order:
# ╠═f4226cfe-ee06-4c72-9615-fc4aedfd045c
# ╟─1f4407d0-9b73-11eb-0e91-cd0de83535aa
