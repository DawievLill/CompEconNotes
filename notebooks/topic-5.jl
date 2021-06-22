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

# ╔═╡ Cell order:
# ╟─eeca0160-d37e-11eb-05f4-bd7469319dad
# ╟─99e0b1b9-d032-41d6-8a1d-742050d5ddeb
# ╟─0ac4cb74-c027-4677-bd45-672f19035ce4
