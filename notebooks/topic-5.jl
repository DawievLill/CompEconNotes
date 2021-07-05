### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ eeca0160-d37e-11eb-05f4-bd7469319dad
begin
	import Pkg
	Pkg.activate(mktempdir())
end

# ╔═╡ dad03d6a-e941-4980-b780-826e417e6311
begin
		Pkg.add([Pkg.PackageSpec(name="Plots")])
		using Plots
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

1. Interpolation or collocation
2. Regression (curve fitting)
3. Local approximations (perturbation methods)

With interpolation we require that an approximation must pass through certain points of the function. The `residual` at each grid point needs to be zero with this class of methods. Finally, in the case of regression we would like to minimise some notion of distance without requiring that the points pass through the function in question. We won't focus too much of our attention on regression problems. In the case of local approximations we can approximate the function and it's derivative at a single point using Taylor series expansions. "

# ╔═╡ 9ffc49e0-056b-4029-87ce-2fee216277a5
md"  ## Basis functions "

# ╔═╡ 31bc8646-bf07-4a22-9030-196af852251c
md" Before we start with the idea of a basis function, we need to understand what a basis is in the context of linear algebra. In general, one of the best investments you can make in your mathematical training is to study linear algebra. One can never know enough linear algebra. For a good description of linear algebra concepts with code, you can visit the [QuantEcon](https://python.quantecon.org/linear_algebra.html) page. We start with the idea of span and the move to linear independence and finally talk about the definition of a basis.

If we are provided a set of vectors $A:=\left\{a_{1}, \ldots, a_{k}\right\}$ in $\mathbb{R}^{n}$ we want to know which new vectors can be created by performing **linear** operations. New vectors that are created in this way are referred to as *linear combinations* of $A$. We have that $y \in \mathbb{R}^{n}$ is a linear combination of $A$ if 

$y=\beta_{1} a_{1}+\cdots+\beta_{k} a_{k} \quad \text{for some scalars} \quad \beta_{1}, \ldots, \beta_{k}$

These $\beta$ values are the coefficients of the linear combination. The set of linear combinations of $A$ is referred to as the **span** of $A$

The following figure shows the span of $A=\left\{a_{1}, a_{2}\right\}$ in $\mathbb{R}^{3}$. It provides a 2 dimensional plane passing through these two points and the origin. 
"

# ╔═╡ 8313e1ee-7c01-4eda-8d4b-888f04de5c90
plotly()

# ╔═╡ b6f096f4-9857-4d66-a62d-8dd818c4d9a2
begin

	# fixed linear function, to generate a plane
	f(x, y) = 0.1x + 0.1y
	
	# lines to vectors
	x_vec = [0 0; 3 3]
	y_vec = [0 0; 3 -3]
	z_vec = [0 0; f(3, 3) f(3, -3)]
	
	# draw the plane
	n = 20
	grid = range(-5, 5, length = n)
	z2 = [ f(grid[row], grid[col]) for row in 1:n, col in 1:n ]
	plot(grid, grid, z2, fill = :ice, st = :surface, alpha = 0.75)
	plot!(x_vec, y_vec, z_vec, colour = [:black :black], linewidth = 5, labels = "", colorbar = false)
	
	
end

# ╔═╡ 7fc20a57-f784-4336-a71f-bb2b0bbdc777
md" For our set of vectors to have a large span, we need linear independence. A collection of vectors is linearly dependent if a strict subset of $A$ has the same span as $A$. The collection of vectors are **linearly independent** if they are not linearly dependent. In other words, no vector is redundant to the span. In our example above, if there were a third vector $a_3$ then the set $\left\{a_{1}, a_{2}, a_{3}\right\}$ would be linearly dependent if $a_3$ lies in the plane. The more formal definition of linear independence is the following,

$\text {If} \quad \beta_{1} a_{1}+\cdots \beta_{k} a_{k}=0 \quad \text{for scalars} \quad \beta_{1}, \ldots, \beta_{k}, \quad \text{then} \quad \beta_{1}=\cdots=\beta_{k}=0$

A basis for the vector space over the real numbers is a linearly independent subset of the vector space that spans that space. In other words, a basis satisfies linear independence and spans the space. Now that we have defined a basis in terms of linear algebra, we need to think about what a **basis function** is. 

Similar to the way in which every vector in a vector space can be represented by a linear combination of basis vectors, we can represent every continuous function in a function space by a linear combination of basis functions. In other words, basis functions are linearly independent functions that span the function space. We are mostly interested in the space of continuous or continuously differentiable functions (our function space of interest).
"


# ╔═╡ 55cfd73a-7015-4936-b229-ddcc488d48de
md" ## Interpolation based methods "

# ╔═╡ 59dffcb5-b978-495e-a25e-22a57ae84882
md" Interpolation is a form of function approximation in which the interpolant and the true function must agree. In other words, they must have same value at a finite number of points. Additional restrictions can also be imposed on the interpolant, such as constraints on monotonicity, convexity and smoothness. In general, interpolation is any method that takes information at a finite set of points and finds a function that satisfies that information.

There are are two primary approaches to interpolation. The approaches are **finite element methods** and **spectral methods**. We will talk about finite element methods first and then move on to spectral methods toward the end of the session. "

# ╔═╡ ef709d7b-107b-43bb-9026-c7af1efacfdc
md" ### Interpolation basics "

# ╔═╡ 5e55cf7d-a076-45e1-8c54-c561be82857d
md" ### Finite element methods "

# ╔═╡ 95125a71-1691-4009-acb4-5db30cd6b69f
md" ### Spectral methods "

# ╔═╡ ef9fb41f-b815-4e77-8913-82ee14fc3d8c
md" This section links interpolation and regression to the concept of a basis function. We will deal with local approximation through Taylor series expansion toward the end of the session. This section draws from Chapter 6 of Judd and Chapter 6 of Miranda and Fackler. 



The space of continuous functions is spanned by monomials, $x^{n}$, $n = 0, 1, 2, \ldots$, and one can use the monomials as a basis for the space of continuous functions. However, we have to ask ourselves whether this is a good basis. Normally a good basis for a vector space also has some orthogonality properties. This leads us into the discussion of constructing orthogonal polynomials as bases for our function space. In order to think about orthogonality one needs the concept of an inner product in the vector space. 

Given that $F$ is a space of continuous real valued functions, let us define a inner product operation on that space. 

$<g, h>=\int_{\mathbf{x}} g(x) h(x) w(x) d x$

where $g, h, w \in F$ and $w$ is a weighting function. The pair $\{F,<.,.>\}$ form an inner-product vetor space."

# ╔═╡ 0a9c73e1-68c1-4c17-bd48-98673b44d332
md" ## Local approximations " 

# ╔═╡ Cell order:
# ╟─eeca0160-d37e-11eb-05f4-bd7469319dad
# ╟─dad03d6a-e941-4980-b780-826e417e6311
# ╟─99e0b1b9-d032-41d6-8a1d-742050d5ddeb
# ╟─0ac4cb74-c027-4677-bd45-672f19035ce4
# ╟─25b6e79e-37dc-4a67-8b56-fbaf9da41713
# ╟─a4f24d72-fb51-47cc-827d-c9f3cb367998
# ╟─9ffc49e0-056b-4029-87ce-2fee216277a5
# ╟─31bc8646-bf07-4a22-9030-196af852251c
# ╠═8313e1ee-7c01-4eda-8d4b-888f04de5c90
# ╠═b6f096f4-9857-4d66-a62d-8dd818c4d9a2
# ╟─7fc20a57-f784-4336-a71f-bb2b0bbdc777
# ╟─55cfd73a-7015-4936-b229-ddcc488d48de
# ╟─59dffcb5-b978-495e-a25e-22a57ae84882
# ╟─ef709d7b-107b-43bb-9026-c7af1efacfdc
# ╟─5e55cf7d-a076-45e1-8c54-c561be82857d
# ╟─95125a71-1691-4009-acb4-5db30cd6b69f
# ╟─ef9fb41f-b815-4e77-8913-82ee14fc3d8c
# ╟─0a9c73e1-68c1-4c17-bd48-98673b44d332
