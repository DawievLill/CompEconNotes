### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 1c7f7f74-7c57-11eb-293a-d1be483a7ca0
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
			Pkg.PackageSpec(name="Images"), 
			Pkg.PackageSpec(name="ImageMagick"), 
			Pkg.PackageSpec(name="PlutoUI"), 
			Pkg.PackageSpec(name="HypertextLiteral"), 
			Pkg.PackageSpec(name="ForwardDiff")
			])

	using Images
	using PlutoUI
	using HypertextLiteral
	using LinearAlgebra
	using ForwardDiff
end

# ╔═╡ 7819e032-7c56-11eb-290b-23dc34edfc58
md"# Functions and Transformations"

# ╔═╡ d88705f0-7c57-11eb-1950-bd54523e4a72
md" This section takes inspiration from a course on computational thinking that is presented at MIT, which can be found [here](https://computationalthinking.mit.edu/Spring21/). Much of what we present here has been taken directly from these notes. We will start with basics on functions and arrays, in order to make sure everyone is one the same page with respect to these fundamental concepts. 

Once these topics have been covered, we move to a really cool way in which you can take derivatives, which is called `autodiff`, which is short for automatic differentiation. With this method we can automatically compute **exact** derivatives (up to floating-point error) given only the function itself.  

This method of taking derivatives is a bit different from the numerical differentiation techniques that we will talk about later in the reading group. It is used widely in machine learning and optimisation and has become increasingly popular over the last couple of years. Finally, we will cover some basic foundational concepts in linear algebra, which will be useful for our next session inverses and solutions of linear systems of equations."

# ╔═╡ 45aed8a2-7c59-11eb-3f69-e701041d6a30
md" ## Functions (in Julia)"

# ╔═╡ 7cfa32e4-7c58-11eb-32c0-5760739f6de4
md"""
Before we get started with automatic differentiation, let us just make sure that everyone is up to speed on the basics of functions. Remember that a function is defined as a relation which assigns to each element in the domain a **single element** in the range. A relation is a set of ordered pairs, $(x, y)$. The set of first coordinates is the domain, the set of second coordinates the range of the relation. Therefore a function is simply a mapping from values in the domain to a single value in the range. 

The first type of functions you learned about were **univariate** functions or real-valued functions functions of a single variable e.g. $f₁(x)=x^2$, $f₂(x)=\sin(x)$, $f₃(x)=x^\alpha$.

In _Julia_, functions can be written in short form, anonymous form, or long form.
"""

# ╔═╡ a3d72150-87e9-11eb-29de-7d2adebbf9a7
f₁(x) = x^2 # Short form

# ╔═╡ cb0ac2ce-87e9-11eb-3c26-ed50a59978be
x -> sin(x) # Anonymous form

# ╔═╡ dfcad438-87e9-11eb-0fd2-218a4806a08c
# Long form (most frequently used)

function f₃(x, α = 3)
	x^α
end

# ╔═╡ 6ad1feb0-87ea-11eb-14fb-2b73c5bacf7d
md" In terms of these univariate functions it is easy to perform automatic differentiation. This method of differentiation is different from symbolic differentiaion that you will encounter in calculus or numerical differentition via first differences."

# ╔═╡ c0cb6d4a-87ec-11eb-348b-e540882173e3
md" There are several packages available to perform automatic differentiation in _Julia_. [Here](https://juliadiff.org/) is a curated list of all actively maintained `autodiff` packages." 

# ╔═╡ b930ded2-87ea-11eb-3f26-3d4e1597615b
md" First, let us illustrate how symbolic differentiation would look in _Julia_ and then turn to automatic differentiation (with numerical methods covered in later sessions)."

# ╔═╡ 22845838-7c5a-11eb-2206-55a258d0d8ee
md" ## Arrays in Julia "

# ╔═╡ b17e4bca-7c5a-11eb-08fa-571bb45b7a3e
md" An array is a rectangular grid that is used for storing data of any type. We can store and retrieve data within this array by using indexing. You might have encountered arrays in mathematics under the name of matrices. An array can be one-dimensional, two-dimensional, three-dimensional and so forth. The dimension gives us an idea of the number of indices that we need to specify. For array objects we also need to know the length of the data in each of the dimensions. "

# ╔═╡ 20c57a22-7c5a-11eb-2b4e-57ef185f5c53
md" It might sound arbitrary at first to focus on things like dimension and length of arrays, but take it from me, you will frequently want to know the properties of the arrays that you are working with."

# ╔═╡ 7ec6d044-7c58-11eb-112d-5d8be6b4288c


# ╔═╡ 7ee10e48-7c58-11eb-2e3b-f5805d4af19c


# ╔═╡ 7efab6f4-7c58-11eb-006e-412089f351a0


# ╔═╡ 7f14d2fa-7c58-11eb-1e88-773a06148b22


# ╔═╡ 7f2f1b06-7c58-11eb-038e-15bd2b4d1dbb


# ╔═╡ 7f495aca-7c58-11eb-0e35-1917e679988c


# ╔═╡ 7f613bb8-7c58-11eb-259f-512e78668dc9


# ╔═╡ 7f7b1022-7c58-11eb-3489-1f1311c1dc44


# ╔═╡ 7f9507f4-7c58-11eb-045d-ad36f4fb184a


# ╔═╡ 7faea394-7c58-11eb-2529-c3881c14b364


# ╔═╡ Cell order:
# ╟─1c7f7f74-7c57-11eb-293a-d1be483a7ca0
# ╟─7819e032-7c56-11eb-290b-23dc34edfc58
# ╟─d88705f0-7c57-11eb-1950-bd54523e4a72
# ╟─45aed8a2-7c59-11eb-3f69-e701041d6a30
# ╟─7cfa32e4-7c58-11eb-32c0-5760739f6de4
# ╠═a3d72150-87e9-11eb-29de-7d2adebbf9a7
# ╠═cb0ac2ce-87e9-11eb-3c26-ed50a59978be
# ╠═dfcad438-87e9-11eb-0fd2-218a4806a08c
# ╟─6ad1feb0-87ea-11eb-14fb-2b73c5bacf7d
# ╟─c0cb6d4a-87ec-11eb-348b-e540882173e3
# ╟─b930ded2-87ea-11eb-3f26-3d4e1597615b
# ╟─22845838-7c5a-11eb-2206-55a258d0d8ee
# ╟─b17e4bca-7c5a-11eb-08fa-571bb45b7a3e
# ╟─20c57a22-7c5a-11eb-2b4e-57ef185f5c53
# ╠═7ec6d044-7c58-11eb-112d-5d8be6b4288c
# ╠═7ee10e48-7c58-11eb-2e3b-f5805d4af19c
# ╠═7efab6f4-7c58-11eb-006e-412089f351a0
# ╠═7f14d2fa-7c58-11eb-1e88-773a06148b22
# ╠═7f2f1b06-7c58-11eb-038e-15bd2b4d1dbb
# ╠═7f495aca-7c58-11eb-0e35-1917e679988c
# ╠═7f613bb8-7c58-11eb-259f-512e78668dc9
# ╠═7f7b1022-7c58-11eb-3489-1f1311c1dc44
# ╠═7f9507f4-7c58-11eb-045d-ad36f4fb184a
# ╠═7faea394-7c58-11eb-2529-c3881c14b364
