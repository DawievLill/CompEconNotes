### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ eeca0160-d37e-11eb-05f4-bd7469319dad
begin
	using ApproxFun
	using BasisMatrices
	using ChebyshevApprox
	using Dierckx
	using Interpolations
	using LinearAlgebra
	using Parameters
	using Plots
	using PlutoUI
	using Polynomials
	using FastGaussQuadrature
	using SpecialFunctions
end

# ╔═╡ 99e0b1b9-d032-41d6-8a1d-742050d5ddeb
html"""
<style>
  main {
    max-width: 900px;
  }
</style>
"""

# ╔═╡ 0a31eff6-d5b0-47c2-8483-042770a77478
md" The following packages were used for the session. Please star 
⭐ the packages on Github to show appreciation. "

# ╔═╡ 0ac4cb74-c027-4677-bd45-672f19035ce4
md" # Function approximation "

# ╔═╡ fba8e77d-9141-4f57-afdd-141bb1311a22
md"""

!!! warning "Reference"
	A good reference for function approximation is the book by [Toby Driscoll](http://tobydriscoll.net/fnc-julia/frontmatter.html). These notes draw heavily form this book.

"""

# ╔═╡ 25b6e79e-37dc-4a67-8b56-fbaf9da41713
md" In many of our economic applications we are often confronted with functions that don't have an analytic representation or are difficult to compute. This will be true for the value functions in our dynamic programming discussion in one of the following sessions.

In order to numerically represent problematic functions on our computer we use approximations, which utilises a simpler function to approximate the more difficult one. Approximations will use data on the function in question. Normally we know some of the function values that correspond to some input points. 

The **goal** is to take our approximation to the data and determine what the function value is on some point that is outside of the finite set of points that we originally evaluated the function on. In economics we do this all the time. We make predictions based on observed data. The big difference with approximation is that we can choose the number of points at which we want to evaluate the function. We produce the data in function approximation. "

# ╔═╡ a4f24d72-fb51-47cc-827d-c9f3cb367998
md" There are three general approaches to approximation (we will only focus on interpolation). 

1. **Interpolation**
2. Regression (curve fitting)
3. Local approximations (perturbation methods)

With interpolation we require that an approximation must pass through certain points of the function. The `residual` at each grid point needs to be zero with this class of methods. 

In the case of regression we would like to minimise some notion of distance without requiring that the points pass through the function in question. We won't focus too much of our attention on regression problems. 

In the case of local approximations we can approximate the function and it's derivative at a single point using Taylor series expansions. "

# ╔═╡ 91fa7272-3e83-45c0-b787-87904d4381fb
md" ## The interpolation problem "

# ╔═╡ 29f3f1a7-c6c1-427a-8781-86ee7dbcd219
md"""

!!! note "The interpolation problem"
	Given $n+1$ distinct points $(t_0,y_0)$, $(t_1,y_1),\ldots,(t_n,y_n)$, with $t_0<t_1<\ldots <t_n$ called **nodes**, the **interpolation problem** is to find a function $\hat f(x)$, called the **interpolant**, such that $\hat f(t_k)=y_k$ for $k=0,\dots,n$. -- from  [FNC (Chapter 5)](http://tobydriscoll.net/fnc-julia/localapprox/interpolation.html)

"""

# ╔═╡ 7e662cb3-8500-4381-ad13-0469f0dbf7eb
md"""

To solve this problem we would naturally suggest polynomials as the best way to approximate the function. Polynomials have nice properties and are easy to work with. However, polynomial interpolation can be a tricky business, as we will soon see. 

"""

# ╔═╡ 2cd62a7a-7573-4511-8c81-488de74ad8f5
md"""

Consider the following example. We are provided the following observations of an unknown function on $[-1, 1]$.

"""

# ╔═╡ 3eaa46e2-5847-4b4d-b48d-fe1142545a8b
md"""
nodes
= $(@bind ν Slider(4:20, show_value=true, default=4))
"""

# ╔═╡ 87a1c372-1e7f-4d00-b897-7739ffdf19e8
begin
	t = range(-1, 1, length = ν + 1)
	y = @. t^2 + t + 0.05 * sin(20 * t) # Broadcast the dot operator with a macro
end

# ╔═╡ 9771216f-9c98-42c7-98a9-ebc264a5ec0c
scatter(t, y, label="data", leg=:top, color = :red)

# ╔═╡ db94347c-df46-43ee-a7c0-0259c11c859f
md"""

We can compute the polynomial interpolant using the `fit` function in Julia. Let us see what this does for our generated data. 

"""

# ╔═╡ ef77b7b9-c2ee-4222-a926-5672678aedeb
begin
	p = Polynomials.fit(t, y, ν)     # interpolating polynomial
	plot!(p, -1, 1, label="interpolant", color = :steelblue, lw = 2)
end

# ╔═╡ 770e4972-6927-4370-8258-25fc101461af
md"""

Shifting around the slider we see that as we increase the number of points, the polynomial interpolant exhibits some strange behaviour near the endpoints. In fact, interpolation by a polynomial at equally spaced nodes is ill-conditioned as the degree of the polynomial grows. 

"""


# ╔═╡ 9ffc49e0-056b-4029-87ce-2fee216277a5
md"  ## Basis (functions) "

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

Similar to the way in which every vector in a vector space can be represented by a linear combination of basis vectors, we can represent every continuous function in a function space by a linear combination of basis functions. In other words, basis functions are linearly independent functions that span the function space. 
"


# ╔═╡ 55cfd73a-7015-4936-b229-ddcc488d48de
md" ## Global methods "

# ╔═╡ 59dffcb5-b978-495e-a25e-22a57ae84882
md" Interpolation is a form of global function approximation in which the interpolant and the true function must agree. In other words, they must have same value at a finite number of points. Additional restrictions can also be imposed on the interpolant, such as constraints on monotonicity, convexity and smoothness. In general, interpolation is any method that takes information at a finite set of points and finds a function that satisfies that information.

There are are two primary approaches to selecting basis functions for interpolation. The approaches are **finite element methods** and **spectral methods**. We will talk about spectral methods first and then move on to finite element methods toward the end of the session. "

# ╔═╡ ef709d7b-107b-43bb-9026-c7af1efacfdc
md" ### Interpolation basics "

# ╔═╡ 2e0a22df-6eab-465b-b23b-a9b14a093290
md" We want to approximate a known function $f(x)$. The interpolant $\hat{f}$ is chosen to be a linear combination of a collection of basis functions $\left\{\phi_{j}(x)\right\}_{j=1}^{n}$ where $n$ is the degree of interpolation. 

Basis functions, as we mentioned before, are linearly independent functions that span the family of functions chosen for the interpolation. In other words, the function space. We are mostly interested in the space of continuous or continuously differentiable functions (our function space of interest).

Thus, for a given family of basis functions $\left\{\phi_{j}(x)\right\}_{j=1}^{n}$ we have 

$f(x) \simeq \hat{f}(x) \equiv \sum_{j=1}^{n} w_{j} \phi_{j}(x)$

Reduced the problem of characterising an infinite dimensional object to the problem of determining the $n$ weights $\{w_{j}\}$. 

The **first step** is to choose an appropriate basis function. Polynomials of increasing order are often used as basis functions.  

The **second step** is designing and interpolation scheme to replicate properties of original function $f$ that one wishes the interpolant to replicate.  

Given that we have $n$ interpolation nodes and $n$ basis functions, computing basis coefficients reduces to solving:

$\sum_{j=1}^{n} w_{j}\phi_{j}\left(x_{i}\right) =f\left(x_{i}\right)=y_{i} \quad \forall i=1,2, \ldots, n$


One can write the problem in matrix notation as follows $\Phi w=y$ where $\Phi_{i j}=\phi_{j}\left(x_{i}\right)$. Written in full,  

$\Phi=\left[\begin{array}{ccc} \phi_{1}\left(x_{1}\right) & \ldots & \phi_{n}\left(x_{1}\right) \\ \vdots & & \\ \phi_{1}\left(x_{n}\right) & \ldots & \phi_{n}\left(x_{n}\right)\end{array}\right], w=\left[\begin{array}{c} w_{1} \\ \vdots \\ w_{n}\end{array}\right], y=\left[\begin{array}{c} y_{1} \\ \vdots \\ y_{n}\end{array}\right]$

The interpolation scheme is well defined if the interpolation nodes and basis function are chosen such that the interpolation matrix is nonsingular. 

Normally, one interpolates to match the value of the original function at selected interpolation nodes. However, we are not limited to point values. We may also use first (and higher order) derivatives at specified points. 

If we have exactly as many nodes as coefficients then this we have a simple root finding problem $\Phi w - y = 0$. If we have more nodes than coefficients then we have a least squares (curve fitting) problem, such that $w = (\Phi^{\prime} \Phi)^{-1} \Phi^{\prime}y$. "



# ╔═╡ 18830fc6-d17d-45f9-8e21-cd42ef74ea4e
md" #### Selection criteria "

# ╔═╡ 36de0424-f52d-40ca-9c42-dfa05277effa
md" From the previous discussion one might be interested in knowing which interpolation nodes and basis functions to choose. Below are some desirable criteria (things that make a good approximation). 

1. The interpolant should approach the target function by increasing the number of nodes and basis functions as well as number of basis once a polynomial basis is chosen. 

2. The basis coefficients should be quick to compute. Works best with diagonal or orthogonal interpolation matrices. 

3. The interpolant should be easy to work with. In other words, the basis function should have the property of being easily evaluated, differentiated and integrated. "

# ╔═╡ e5910bba-3703-4eb0-9aae-55554193de20
md" #### Note on node placement "

# ╔═╡ 095e025d-0344-4e2c-97d9-fa1b230b2469
md" One option for the interpolation nodes is to select evenly spaced nodes. This yields the following selection for nodes,

$x_{i}=a+\frac{i-1}{n-1}(b-a), \quad \forall i=1,2, \ldots, n$

Evenly spaced nodes are not always a good choice, even if your function is smooth. Theory and empirical evidence suggests that Chebyshev nodes are better to use than evenly spaced nodes. These are the roots of the Chebyshev polynomial (which we will discuss soon). "


# ╔═╡ 75ba3ecf-6dff-482d-99f3-9650b6597540
md" Chebyshev nodes are defined on the interval $[-1,1]$ as 

$x_{i}=\cos \left(\frac{2 k-1}{2 n} \pi\right), k=1, \ldots, n$

These nodes can be mapped to any general interval $[a, b]$ by an appropriate transformation. Below is a graph that shows the Chebyshev nodes on the interval $[-1, 1]$."

# ╔═╡ be3e84ef-842f-43ed-833b-7faee299c904
md"""
nodes = $(@bind m Slider(3:100, show_value=true, default=3))
"""

# ╔═╡ 41cee516-36a5-4af3-b157-46fcd1ff6aa5
begin
	plotly()
	gcnodes = gausschebyshev(m)
	scatter(gcnodes, ones(m), ylims=(0.9,1.1),m=(:red),label = "Chebyshev nodes",size=(600,400))
end

# ╔═╡ 93a84049-7e73-41c4-afd8-0359656c6399
md" One can see that the nodes are not evenly scattered. They are spaced further apart in the center of the interval and then more closely grouped toward the end points of the interval. Next we need to think about what type of basis we can use for the problem. 

General advice for node placement is that we want more nodes where the function changes rapidly or has kinks. Usually these points are close to the endpoints of the interval. Now we move to the selection of basis functions. "

# ╔═╡ 95125a71-1691-4009-acb4-5db30cd6b69f
md" ### Spectral methods "

# ╔═╡ 507628c6-2c13-4ff5-a73f-c6f30c8e5a8c
md" Spectral methods are differentiated from finite element methods in that basis functions are **nonzero over the entire domain** of the function. 

Spectral methods basically refer to the usage of a polynomial basis. Polynomials are generally nonzero. The motivation for using polynomials originates form the famous **Weierstrass Theorem**. This theorem states that there exists a polynomial that approximates any continuous function over a compact (closed and bounded) domain arbitrarily well. This theorem does not tell us which type of polynomial to use. This is something that we need to figure out for ourselves. "

# ╔═╡ b1f875af-4046-445e-a2fb-afca07f992ef
md" There are many polynomial families that are well suited to approximate continuous and differentiable functions. However we will focus mostly on monomials and orthogonal polynomials. In the class of orthogonal polynomials most of our attention will be placed on Chebyshev polynomials. "

# ╔═╡ 8918fde2-07d7-4f76-8c5b-03fa04650058
md" #### Monomials "

# ╔═╡ c4e41e78-e554-47b3-822f-3516f615bc2c
md" The most logical choice for a basis is the monomial basis. The monomials $x^{i}$, $i = 1, 2, \ldots$ build a basis for the space of continuous functions over the interval $[a, b]$. Each member of this function space can be represented by 

$\hat{f}(x)=\sum_{j=1}^{\infty} w_{j} x^{j-1}$

Common to use a linear combination of the first $n$ members of this basis to approximate continuous functions, 

$f(x) \approx \hat{f}(x)=w_{1}+w_{2} x+w_{3} x^{2}+\ldots+w_{n} x^{n-1}$

The figure below shows that the monomials basis components become more similar to each other as we increase the degree of the basis. In other words, the additional information is marginally less valuable."

# ╔═╡ 1d8391b6-c28b-4707-9a27-b6152b14be49
function plot_function(basis_function, x, n)
	
	for i = 1:n-1
		
		f_data = basis_function(x, i)
		
		if i == 1
			plot(x, f_data, linewidth = 1.5, xlabel = "x", ylabel = "Basis functions", label = "", 
			tickfontsize = 10, guidefontsize = 10);
		else
			plot!(x, f_data, linewidth = 1.5, label = "");
		end
	end
	
	f_data = basis_function(x, n)
	plot!(x, f_data, linewidth = 1.5, label = "")
end	# Read through this again to see how the plotting scheme works here. Struggling to make sense of this. 

# ╔═╡ 606948ed-57bd-4f25-8b8a-99133c405bc8
x_range = -1:.01:1;

# ╔═╡ c232b1db-ab1c-4c29-b0ee-8fc719b9d2bc
monomials(x, n) = x.^n

# ╔═╡ 79b677e8-98a7-446e-babb-c9f8ccdec52b
md"""
Polynomial degree = $(@bind τ Slider(2:20, show_value=true, default=2))
"""

# ╔═╡ e45843c7-a0a9-40fd-abff-d150fc1efb67
plot_function(monomials, x_range, τ)

# ╔═╡ 4e21ef20-c5ae-46c7-a512-efb47c48a244
md" This idea is also captured in the basis matrix for monomials. The matrix looks as follows, 

$\left[\begin{array}{cccc} 1 & x_{1} & \ldots & x_{1}^{n} \\ & \vdots & \ddots & \vdots \\ 1 & x_{n} & \ldots & x_{n}^{n}\end{array}\right]\left[\begin{array}{c} w_{1} \\ \vdots \\ w_{n}\end{array}\right]=\left[\begin{array}{c} y_{1} \\ \vdots \\ y_{n}\end{array}\right]$

This matrix is called a Vandermonde matrix. Theorethically this system is well defined, but in practice this matrix is known for being ill-conditioned, especially for high-degree polynomials. The columns of the matrix are nearly linearly dependent as the degree increases. Let us look at an example where we approximate a function with a monomial basis on an evenly spaced grid."

# ╔═╡ 46bf582a-ba9a-491f-96f9-fa1ec38698a9
function monomial(f, n, lb, ub)
	
    # This function solves Φw = y → w = Φ\y
    # Φ = matrix of monomial basis functions evaluated on the grid
    
	coll_points = range(lb, ub, length = n)                      # evenly spaced collocation nodes
    y_values = f(coll_points)                                    # function values on the grid
    basis_functions = [coll_points.^degree for degree = 0:n-1]   # vector of basis functions
    basis_matrix = hcat(basis_functions...)                      # basis matrix
    coefficients = basis_matrix\y_values                         # w = Φ\y
	
    return coefficients
end

# ╔═╡ 53bee8bb-5469-4f99-86f1-a6e9ce0f44f7
function f_approx(coefficients, points)
	
	n = length(coefficients) - 1
	basis_functions = [coefficients[degree + 1] * points.^degree for degree = 0:n]  # evaluate basis functions 
	basis_matrix = hcat(basis_functions...)										   # transform into matrix
	function_values = sum(basis_matrix, dims = 2)								   # sum up into function value
	
	return function_values
end

# ╔═╡ 75885058-267f-4e7e-9cb2-bd5e1f2c87c2
md"""
α = $(@bind α Slider(4:8, show_value=true, default=4))
"""

# ╔═╡ 456bcf5c-25a2-453d-80ac-3b103a671d93
g(x) = sin.(x)

# ╔═╡ 3606489f-e45d-4f06-ab72-2c31def7ed10
plot_points = 0:.01:2pi;

# ╔═╡ 93d60563-a66e-4a8f-a11b-a81564cd3941
coefficients_α = monomial(g, α, 0, 2pi);

# ╔═╡ b089c56d-9f9e-48c8-bf58-611b58bdcd0c
f_values_α = f_approx(coefficients_α, plot_points);

# ╔═╡ 2ca3cba6-90db-4d4b-a225-7ed25093fedd
begin
	plot(g, plot_points, label = "sin(x)", line = 3.5)
	plot!(plot_points, f_values_α, label = "Approximation", line = 2, linestyle = :dash)
end

# ╔═╡ d8b7067c-a4be-4d87-bd79-adb721c61935
md" Like we said before, we will almost never use monomial basis functions. Try this same method on the Runge function $f(x) = 1 / (1 + 25x^2)$ to see why this is the case. "

# ╔═╡ 3bd260fc-de99-4cf9-9453-828df9178b64
md"""
β = $(@bind β Slider(3:20, show_value=true, default=3))
"""

# ╔═╡ ec27fb6b-1cbb-46c8-9da6-bf128dd5feae
runge(x) = 1 ./ (1 .+ 25x.^2)	

# ╔═╡ ebfa42bb-0c34-40f1-82ab-dd40d8954eb9
runge_plot_points = -1:.01:1;

# ╔═╡ a2ef0e9f-8dc9-4892-acbb-3b7b13681c93
runge_coefficients_β = monomial(runge, β, -1, 1);

# ╔═╡ b1f8acd8-3de8-4512-91e1-47b45780c0ea
runge_values_β = f_approx(runge_coefficients_β, runge_plot_points);

# ╔═╡ 38a7ad81-944a-438a-ba07-75c24daffdd1
begin
	plot(runge, runge_plot_points, label = "Runge Function", line = 3)
	plot!(runge_plot_points, runge_values_β, label = "Approximation", line = 2, linestyle = :dash)
end

# ╔═╡ 113e9d8d-17cc-413e-95ca-0c53dda1f48b
md" #### Orthogonal polynomials "

# ╔═╡ ef9fb41f-b815-4e77-8913-82ee14fc3d8c
md" 
The space of continuous functions is spanned by monomials, $x^{n}$, $n = 0, 1, 2, \ldots$, and one can use the monomials as a basis for the space of continuous functions. However, we have to ask ourselves whether this is a good basis. Normally a good basis for a vector space also has some orthogonality properties. This leads us into the discussion of constructing orthogonal polynomials as bases for our function space. Let us introduce some definitions in order to facilitate the discussion. 

**Definition (weighting function)** A weighting function $\omega(x)$ on an interval $[a, b]$ is a function that is positive almost everywhere on $[a, b]$ and has finite integral on $[a, b]$.

An example of a weighting function would be $\omega(x) = (1 - x)^{-1/2}$ over the interval $[-1, 1]$. This function is positive over the entire interval and integrates to a value of $\pi$ over the interval.

In order to think about orthogonality one needs the concept of an inner product in the vector space.

**Definition (inner product)** Given two functions $f_1(x)$ and $f_2(x)$ defined on $[a, b]$, the inner product with respect to the weighting function is given by

$\left\langle f_{1}, f_{2}\right\rangle=\int_{a}^{b} f_{1}(x) f_{2}(x) \omega(x) \mathrm{d} x$

As an example, assume that $f_1 = 1$ and $f_2  = x$ and $\omega = (1 - x)^{-1/2}$. The inner product over the interval $[-1, 1]$ is 

$\left\langle f_{1}, f_{2}\right\rangle=\int_{-1}^{1} \frac{x}{\sqrt{1-x}} \mathrm{~d} x=-\left.\sqrt{1-x^{2}}\right|_{-1} ^{1}=0$

In this case we actually have the orthogonality property with respect to the weighting function. 

**Definition (orthogonal polynomial)** Family of polynomials $\left\{P_{n}(x)\right\}$ is mutually orthogonal with respect to the $\omega(x)$ iff

$\left\langle P_{i}, P_{j}\right\rangle=0 \quad \text{for} \quad i \neq j$

Below is a table of the most common families of orthogonal polynomials. 

$$\begin{aligned}
&\text { Orthogonal polynomials (definitions) }\\
&\begin{array}{lccl}
\hline \hline \text { Family } & \omega(x) & {[a ; b]} & \text { Definition } \\
\hline \text { Legendre } & 1 & {[-1 ; 1]} & P_{n}(x)=\frac{(-1)^{n}}{2^{n} n !} \frac{\mathrm{d}^{n}}{\mathrm{~d} x^{n}}\left(1-x^{2}\right)^{n} \\
\text { Chebychev } & \left(1-x^{2}\right)^{-1 / 2} & {[-1 ; 1]} & T_{n}(x)=\cos \left(n \cos ^{-1}(x)\right) \\
\text { Laguerre } & \exp (-x) & {[0, \infty)} & L_{n}(x)=\frac{\exp (x)}{n !} \frac{\mathrm{d}^{n}}{\mathrm{~d} x^{n}}\left(x^{n} \exp (-x)\right) \\
\text { Hermite } & \exp \left(-x^{2}\right) & (-\infty, \infty) & H_{n}(x)=(-1)^{n} \exp \left(x^{2}\right) \frac{\mathrm{d}^{n}}{\mathrm{~d} x^{n}} \exp \left(-x^{2}\right) \\
\hline \hline
\end{array}
\end{aligned}$$

"

# ╔═╡ 1d3fa479-5844-49ae-bb3c-5a12954cdf60
md" ##### Chebyshev polynomials "

# ╔═╡ 37f7640f-13ab-4d39-a9f1-7a258fe62bf6
md" Orthogonal polynomials span the polynomial space. Below is a figure that provides different degrees of the Chebyshev polynomial over the interval $[-1, 1]$. These polynomials are an excellent basis for constructing polynomials that interpolate function values at the Chebyshev nodes (which were mentioned before).

Chebychev basis polynomials in combination with Chebychev interpolation nodes yields an extremely well-conditioned interpolation equation that can be accurately and effciently solved, even with high
degree approximants. "

# ╔═╡ cb81c6b1-a8ec-473d-bc27-1bdbb46b5a35
md"""
Polynomial degree = $(@bind δ Slider(2:10, show_value=true, default=2))
"""

# ╔═╡ 97b9730f-fffe-4e8d-86e7-553cdd8dcf5b
function cheb_polys(x, n)
    if n == 0
        return x./x                 # T_0(x) = 1
    elseif n == 1
        return x                    # T_1(x) = x
    else
        cheb_recursion(x, n) =
            2x.*cheb_polys.(x, n - 1) .- cheb_polys.(x, n - 2)
        return cheb_recursion(x, n) # T_n(x) = 2xT_{n-1}(x) - T_{n-2}(x)
    end
end;

# ╔═╡ f371ee9e-c72e-4dc0-90d2-b99345edd4d0
plot_function(cheb_polys, x_range, δ)

# ╔═╡ e46bc8e3-10d7-4637-8d68-9b9ad5f9bf41
md" These polynomials are easily calculated through recursion. 

$\begin{align} T_0(x) & = 1 \\
  T_1(x) & = x \\
  T_{n+1}(x) & = 2xT_n(x) - T_{n-1}(x) \end{align}$

The explicit formulation as stated before is $T_{n}(x)=\cos \left(n \cos ^{-1}(x)\right)$. "

# ╔═╡ da986aee-d913-4c23-b214-dd4ae0ede309
function T(x,n)
    @assert (x >= -1) & (x <= 1)
    if n == 0
        return 1.0
    elseif n==1
        return x
    else
        2*x*T(x,n-1) - T(x,n-2)
    end
end

# ╔═╡ d27de7be-46b7-4801-b3c3-920ae3c00efe
T2(x,n) = cos(n* acos(x))

# ╔═╡ 9b1b2070-bda5-43d9-be31-b83829fcd98d
T2(0.5,3) == T(0.5,3) # Compare the two methods. 

# ╔═╡ 59214606-5ad7-4905-b029-43f4271fe5b5
md" Since these polynomials span the entire polynomial vector space, the $\Phi$ matrix has full rank and is invertible. There are many reasons, both theorethical and practical to use these polynomials. There is a famous quote by John Boyd on the three moral principles for the selection of Chebyshev polynomials. 

1. When in doubt, use Chebyshev polynomials unless the solution is spatially periodic, in which case an ordinary Fourier series is better
2. Unless you are sure another set of basis functions is better, use Chebyshev polynomials
3. Unless you are really, really sure another set of basis functions is better use Chebyshev polynomials

The main message is that for periodic functions, use Fourier series. For everything else, use Chebyshev polynomials.  "

# ╔═╡ a96fe789-13d8-42f3-9c04-07e56637949c
md" #### ApproxFun.jl "

# ╔═╡ f7dba6aa-e1fd-451c-9ab7-e76a483f6959
md" One of the packages that we can use for function approximation and manipulation is ApproxFun.jl. This package allows for more than simple function approximation functionality. You can do algebra, differentiate and integrate functions, solve ODES and PDEs, represent periodic functions, etc. This is an amazing package and we won't have time cover most of3 it in this session, but I highly recommend you go to the website and take a look. It draws inspiration from the Matlab package chebfun, which is quite famous. Nick Hale at the applied mathematics department worked on chebfun when he was at Oxford, so you can ask him some questions on the package. "

# ╔═╡ 2fa769ef-29d7-440a-89fb-c24980e609ca
z = Fun(runge,-1..1)

# ╔═╡ f0786957-9371-4569-9631-e173b068530e
space(z)

# ╔═╡ a6837eae-68ac-420b-b47e-a05a32e355a5
begin
	plot(z; label = "Runge Approx", line = 3, linestyle = :dash)
	plot!(runge; label = "Runge Function", line = 2, alpha = 0.7)
end

# ╔═╡ 70edc423-3403-4652-b7a6-e90720543918
ramp(x) = (x + abs(x)) / 2

# ╔═╡ 11b12320-f60e-4a6f-a2fb-b85a573abdd2
q = Fun(ramp, -1..1)

# ╔═╡ 7387b11f-873b-4af1-b909-bf3051d16e6e
begin
	plot(q; label = "Ramp Approx", line = 3, linestyle = :dash)
	plot!(ramp; label = "Ramp Function", line = 2, alpha = 0.7)
end

# ╔═╡ 5e55cf7d-a076-45e1-8c54-c561be82857d
md" ### Finite element methods "

# ╔═╡ e5ab6234-6221-4ef2-b2b1-05f7e8751928
md" We can also map data values to functions via finite element methods. With these methods, the basis functions are nonzero over **subintervals** of the domain of the function. Finite element methods include different versions of piecewise polynomials over segments of the domain. Our discussion will focus on different types of splines. "

# ╔═╡ d74faeb0-2338-4f9d-9cae-a240afce8c83
md" #### Piecewise linear interpolation "

# ╔═╡ fc1ca90d-acbe-4a59-984f-4f9afdc18c42
md" We look at the same interpolation example from above and this time consider different finite element methods. Piecewise linear interpolation is simply a game of connect-the-dots. That is, the data points are joined pairwise by line segments."

# ╔═╡ 01979f83-69f4-470a-9f47-8d89663b9b9e
md" 

Given certain nodes $t_0 < t_1 < \cdots < t_n$, the piecewise linear interpolant $p(x)$ is given by, 

$$p(x) = y_k + \frac{y_{k+1}-y_k}{t_{k+1}-t_k}(x-t_k) \quad \text{ for } x\in[t_k,t_{k+1}]$$

Here we have that in each interval $[t_k, t_{k+1}]$ the interpolant is a linear function that passes through both the points $(t_k, y_k)$ and $(t_{k+1}, y_{k+1})$. 

We have already discussed the idea of a basis of the set of functions. We can choose the interpolant from among the linear combination of a preselected finite set of functions. In the current setting we use 

$$\begin{split}  H_k(x) =
  \begin{cases}
    \dfrac{x-t_{k-1}}{t_k-t_{k-1}} & \text{if $x\in[t_{k-1},t_k]$},\\[2.5ex]
    \dfrac{t_{k+1}-x}{t_{k+1}-t_{k}} & \text{if $x\in[t_{k},t_{k+1}]$},\\[2.5ex]
    0 & \text{otherwise}.
  \end{cases} \qquad  \end{split}$$

The functions $H_0,\ldots,H_n$ are called **hat functions**. Each of these hat functions is globally continuous and linear inside the interval $[t_k, t_{k+1}]$. 

Any linear combination of these hat functions also has the property of being continuous. 

Any hat function is expressible as a unique linear combination of hat functions, 

$$\sum_{k=0}^n c_k H_k(x)$$

Hat functions therefore form a basis for the set of functions that are continuous and piecewise linear. 

"

# ╔═╡ 62a9d6b6-1950-4866-9e2f-86fde888060e
md"""
points = $(@bind ν_1 Slider(4:20, show_value=true, default=4))
"""

# ╔═╡ 6291861f-e52a-4c74-a0bf-fed791d54ef1
begin
	t_1 = range(-1, 1, length = ν_1 + 1)
	y_1 = @. t_1^2 + t_1 + 0.05 * sin(20 * t_1) # Broadcast the dot operator with a macro
end

# ╔═╡ 65a2c53c-3a15-4f48-8cb7-7eb3b3c1b6fb
scatter(t_1, y_1, label="data", leg=:top, color = :red)

# ╔═╡ 3cf71a1d-0566-4f10-9e71-12f4b9b92eb8
md" We can easily try and connect the dots with linear interpolation from the `Interpolations.jl` package."

# ╔═╡ 9e76aba3-ba12-4212-ad51-2e0b0dd0dee7
begin
	itp_linear = LinearInterpolation(t_1, y_1)
	f_linear(t_1) = itp_linear(t_1)
end

# ╔═╡ 6c4c01c4-78a0-4a7b-9268-743f8fcf763b
begin
	scatter(t_1, y_1, markersize=3,label="Data points")
	plot!(f_linear, t_1, w=1.5,label="Piecewise linear", leg = :top)
end

# ╔═╡ 02edea51-88dd-4ead-af57-9cf844e70d82
md" #### Cubic interpolation "

# ╔═╡ 19962772-9ddf-4c16-bb2a-bee7abd041e1
md" Another option here is to use piecewise cubic interpolation using the `Dierckx.jl` package"

# ╔═╡ 41a2e257-699a-4c1f-b6de-3441aa0a3cfc
begin
	p_1 = Spline1D(t_1,y_1)
	plot!(x->p_1(x),-1,1,label="Piecewise cubic")
end

# ╔═╡ 0a9c73e1-68c1-4c17-bd48-98673b44d332
md" ## Collocation " 

# ╔═╡ 5dad0cee-8239-43df-b1df-bf47b2d7ec8b


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ApproxFun = "28f2ccd6-bb30-5033-b560-165f7b14dc2f"
BasisMatrices = "08854c51-b66b-5062-a90d-8e7ae4547a49"
ChebyshevApprox = "17a596ad-87cd-578c-9b6d-01108c31dc04"
Dierckx = "39dd38d3-220a-591b-8e3c-4c3a8c710a94"
FastGaussQuadrature = "442a2c76-b920-505d-bb47-c5924d526838"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Polynomials = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[compat]
ApproxFun = "~0.12.6"
BasisMatrices = "~0.7.0"
ChebyshevApprox = "~0.1.12"
Dierckx = "~0.5.2"
FastGaussQuadrature = "~0.4.5"
Interpolations = "~0.13.6"
Parameters = "~0.12.2"
Plots = "~1.20.1"
PlutoUI = "~0.7.9"
Polynomials = "~1.2.1"
SpecialFunctions = "~0.10.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[ApproxFun]]
deps = ["AbstractFFTs", "ApproxFunBase", "ApproxFunFourier", "ApproxFunOrthogonalPolynomials", "ApproxFunSingularities", "Calculus", "DomainSets", "DualNumbers", "FFTW", "FastTransforms", "LinearAlgebra", "RecipesBase", "Reexport", "SpecialFunctions"]
git-tree-sha1 = "cada9d687196534b6352543cd6a5004f60a5d646"
uuid = "28f2ccd6-bb30-5033-b560-165f7b14dc2f"
version = "0.12.6"

[[ApproxFunBase]]
deps = ["AbstractFFTs", "BandedMatrices", "BlockArrays", "BlockBandedMatrices", "Calculus", "DSP", "DomainSets", "DualNumbers", "FFTW", "FastGaussQuadrature", "FillArrays", "InfiniteArrays", "IntervalSets", "LazyArrays", "LinearAlgebra", "LowRankApprox", "SparseArrays", "SpecialFunctions", "StaticArrays", "Statistics", "Test", "ToeplitzMatrices"]
git-tree-sha1 = "bab00664565c82f38f562372951ca9b131072b5b"
uuid = "fbd15aa5-315a-5a7d-a8a4-24992e37be05"
version = "0.3.14"

[[ApproxFunFourier]]
deps = ["AbstractFFTs", "ApproxFunBase", "DomainSets", "FFTW", "FastTransforms", "InfiniteArrays", "IntervalSets", "LinearAlgebra", "Reexport"]
git-tree-sha1 = "3209c243eff4ced9d4b422512d23d0a22b2b4cdd"
uuid = "59844689-9c9d-51bf-9583-5b794ec66d30"
version = "0.2.15"

[[ApproxFunOrthogonalPolynomials]]
deps = ["AbstractFFTs", "ApproxFunBase", "BandedMatrices", "BlockArrays", "BlockBandedMatrices", "DomainSets", "FFTW", "FastGaussQuadrature", "FastTransforms", "FillArrays", "IntervalSets", "LinearAlgebra", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "a7329c5310e4c70e93516b96756c6cf5beafe905"
uuid = "b70543e2-c0d9-56b8-a290-0d4d6d4de211"
version = "0.3.10"

[[ApproxFunSingularities]]
deps = ["ApproxFunBase", "ApproxFunOrthogonalPolynomials", "DomainSets", "IntervalSets", "LinearAlgebra", "Reexport", "Statistics"]
git-tree-sha1 = "77831f1f073aeebc7ea56bace5e3e02e07828e23"
uuid = "f8fcb915-6b99-5be2-b79a-d6dbef8e6e7e"
version = "0.2.5"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "1ee88c4c76caa995a885dc2f22a5d548dfbbc0ba"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.2.2"

[[ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "40debc9f72d0511e12d817c7ca06a721b6423ba3"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.17"

[[ArrayLayouts]]
deps = ["Compat", "FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "8f6af27c33b766f19fa6cfe46e629775cda81f88"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.4.11"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[BandedMatrices]]
deps = ["ArrayLayouts", "Compat", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "8c83ee44dc9835263760ad4e77ed4eed4b3490c1"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.15.25"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BasisMatrices]]
deps = ["Combinatorics", "LinearAlgebra", "QuantEcon", "SparseArrays", "Statistics"]
git-tree-sha1 = "d0206f53807fa78720b56f04bb38ef111ec8c7be"
uuid = "08854c51-b66b-5062-a90d-8e7ae4547a49"
version = "0.7.0"

[[BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[BinaryProvider]]
deps = ["Libdl", "Logging", "SHA"]
git-tree-sha1 = "ecdec412a9abc8db54c0efc5548c64dfce072058"
uuid = "b99e7846-7c00-51b0-8f62-c81ae34c0232"
version = "0.5.10"

[[BlockArrays]]
deps = ["ArrayLayouts", "Compat", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "824b1094a47d7da81f9ff77cb56c3341f2f92097"
uuid = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
version = "0.12.14"

[[BlockBandedMatrices]]
deps = ["ArrayLayouts", "BandedMatrices", "BlockArrays", "Distributed", "FillArrays", "LinearAlgebra", "MatrixFactorizations", "Pkg", "SharedArrays", "SparseArrays", "Statistics"]
git-tree-sha1 = "50374b927844af8c21d31db10b833af39493eaca"
uuid = "ffab5731-97b5-5995-9138-79e8c1846df0"
version = "0.9.5"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "80ca332f6dcb2508adba68f22f551adb2d00a624"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.3"

[[ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[ChebyshevApprox]]
deps = ["LinearAlgebra", "ThreadPools"]
git-tree-sha1 = "19d88cd5fbc46bf2893873120763a31c6d351a86"
uuid = "17a596ad-87cd-578c-9b6d-01108c31dc04"
version = "0.1.12"

[[CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "78bee250c6826e1cf805a88b7f1e86025275d208"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.46.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DSP]]
deps = ["FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "2a63cb5fc0e8c1f0f139475ef94228c7441dc7d0"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.6.10"

[[DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Dierckx]]
deps = ["Dierckx_jll"]
git-tree-sha1 = "633c119fcfddf61fb4c75d77ce3ebab552a44723"
uuid = "39dd38d3-220a-591b-8e3c-4c3a8c710a94"
version = "0.5.2"

[[Dierckx_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6596b96fe1caff3db36415eeb6e9d3b50bfe40ee"
uuid = "cd4c43a9-7502-52ba-aa6d-59fb2a88580b"
version = "0.1.0+0"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "992a23afdb109d0d2f8802a30cf5ae4b1fe7ea68"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.11.1"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "StaticArrays", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "501c11d708917ca09ce357bed163dbaf0f30229f"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.23.12"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[DomainSets]]
deps = ["IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics", "Test"]
git-tree-sha1 = "372f681128afe9e728ff00dca2df6b0795e51fbe"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.4.2"

[[Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "03b753748fd193a7f2730c02d880da27c5a24508"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.6.0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "ccd479984c7838684b3ac204b716c89955c76623"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+0"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions"]
git-tree-sha1 = "6ea5f7b4aecce0e3a14ca1da03f62f86148c8fa3"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "0.4.5"

[[FastTransforms]]
deps = ["AbstractFFTs", "ArrayLayouts", "BinaryProvider", "DSP", "FFTW", "FastGaussQuadrature", "FastTransforms_jll", "FillArrays", "Libdl", "LinearAlgebra", "Reexport", "SpecialFunctions", "Test", "ToeplitzMatrices"]
git-tree-sha1 = "cc85cb0f94ba4eee6b0c635cb0c48d54b8bed428"
uuid = "057dd010-8810-581a-b7be-e3fc3b93f78c"
version = "0.12.6"

[[FastTransforms_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "FFTW_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl", "MPFR_jll", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "21dc8ae16745cafc5708b8615e57443f0bed62fd"
uuid = "34b6f7d7-08f9-5794-9e10-3819e4c7e49a"
version = "0.5.4+1"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "502b3de6039d5b78c76118423858d981349f3823"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.9.7"

[[FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "e3af8444c9916abed11f4357c2f59b6801e5b376"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.13.1"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "187198a4ed8ccd7b5d99c41b69c679269ea2b2d4"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.32"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[GMP_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "781609d7-10c4-51f6-84f2-b8444358ff6d"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "182da592436e287758ded5be6e32c406de3a2e47"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.58.1"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "2d908286d120c584abbe7621756c341707096ba4"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.66.2+0"

[[GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "fb28b5dc239d0174d7297310ef7b84a11804dfab"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.0.1"

[[GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "a7a97895780dab1085a97769316aa348830dc991"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.3"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[InfiniteArrays]]
deps = ["DSP", "FillArrays", "LazyArrays", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "aa50bca542c26505efb399290eb5f83e995517bc"
uuid = "4858937d-0d70-526a-a4dd-2d5cb5dd786c"
version = "0.8.2"

[[Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "d19f9edd8c34760dca2de2b503f969d8700ed288"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.4"

[[IntegerMathUtils]]
git-tree-sha1 = "f366daebdfb079fd1fe4e3d560f99a0c892e15bc"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.0"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b7bc05649af456efc75d178846f47006c2c4c3c7"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.6"

[[IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "bcf640979ee55b652f3b01650444eb7bbe3ea837"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.4"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "f3c7f871d642d244e7a27e3fb81e8441e13230d8"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.8.0"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ad927676766e6529a2d5152f12040620447c0c9b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "14.0.4+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "1a43be956d433b5d0321197150c2f94e16c0aaa0"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.16"

[[LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "MatrixFactorizations", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "1e9f6f50e6b39b2cabb18d5f0fafdd45d9c2a28f"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "0.18.1"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LightGraphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "432428df5f360964040ed60418dd5601ecd240b6"
uuid = "093fc24a-ae57-5d10-9952-331d41423f4d"
version = "1.3.5"

[[LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LowRankApprox]]
deps = ["FFTW", "FillArrays", "LinearAlgebra", "Nullables", "Random", "SparseArrays", "Test"]
git-tree-sha1 = "35e31d7e505492cb70b44eacf94e89adf6ef79f6"
uuid = "898213cb-b102-5a47-900c-97e73b919f73"
version = "0.4.3"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "41d162ae9c868218b1f3fe78cba878aa348c2d26"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.1.0+0"

[[MPFR_jll]]
deps = ["Artifacts", "GMP_jll", "Libdl"]
uuid = "3a97d323-0669-5f0c-9066-3539efd106a3"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "JSON", "LinearAlgebra", "MutableArithmetics", "OrderedCollections", "Printf", "SparseArrays", "Test", "Unicode"]
git-tree-sha1 = "23c99cadd752cc0b70d4c74c969a679948b1bb6a"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.2.0"

[[MathProgBase]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9abbe463a1e9fc507f12a69e7f29346c2cdc472c"
uuid = "fdba3010-5040-5b88-9595-932c9decdf73"
version = "0.7.8"

[[MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Random"]
git-tree-sha1 = "292e5f9f0761f3511edfb5420b4feadd9ba165b0"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "0.6.1"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "ae6676d5f576ccd21b6789c2cbe2ba24fcc8075d"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.5"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "29714d0a7a8083bba8427a4fbfb00a540c681ce7"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.3"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "4e675d6e9ec02061800d6cfb695812becbd03cdf"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.4"

[[NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[NLopt]]
deps = ["MathOptInterface", "MathProgBase", "NLopt_jll"]
git-tree-sha1 = "5a7e32c569200a8a03c3d55d286254b0321cd262"
uuid = "76087f3c-5699-56af-9a33-bf431cd00edd"
version = "0.6.5"

[[NLopt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9b1f15a08f9d00cdb2761dcfa6f453f5d0d6f973"
uuid = "079eb43e-fd8e-5478-9966-2cf3e3edb778"
version = "2.7.1+0"

[[NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Nullables]]
git-tree-sha1 = "8f87854cc8f3685a60689d8edecaa29d2251979b"
uuid = "4d1e1d77-625e-5b40-9113-a560ec7a8ecd"
version = "1.0.0"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "1ea784113a6aa054c5ebd95945fa5e52c2f378e7"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.7"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e60321e3f2616584ff98f0a4f18d98ae6f89bbb3"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.17+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "ad8de074ed5dad13e87d76c467a82e5eff9c693a"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.2"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse", "Test"]
git-tree-sha1 = "95a4038d1011dfdbde7cecd2ad0ac411e53ab1bc"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.10.1"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "2276ac65f1e236e0a6ea70baff3f62ad4c625345"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.2"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "3d5bf43e3e8b412656404ed9466f1dcbf7c50269"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.0"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "9888e59493658e476d3073f1ce24348bdc086660"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.0"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "8365fa7758e2e8e4443ce866d6106d8ecbb4474e"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.20.1"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "OffsetArrays", "RecipesBase"]
git-tree-sha1 = "0b15f3597b01eb76764dd03c3c23d6679a3c32c8"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "1.2.1"

[[PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "311a2aa90a64076ea0fac2ad7492e914e6feeb81"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.3"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "3c009334f45dfd546a16a57960a821a1a023d241"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.5.0"

[[QuantEcon]]
deps = ["DSP", "DataStructures", "Distributions", "FFTW", "LightGraphs", "LinearAlgebra", "Markdown", "NLopt", "Optim", "Pkg", "Primes", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "Test"]
git-tree-sha1 = "4e2dc3044303aa2cbf6e321cb9af3982f6774e6a"
uuid = "fcd29c91-0bd7-5a09-975d-7ac3f643a60c"
version = "0.16.2"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "1f27772b89958deed68d2709e5f08a5e5f59a5af"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.3.7"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["OpenSpecFun_jll"]
git-tree-sha1 = "d8d8b8a9f4119829410ecd706da4cc8594a1e020"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "0.10.3"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "7f5a513baec6f122401abfc8e9c074fdac54f6c1"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.1"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "da4cf579416c81994afd6322365d00916c79b8ae"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "0.12.5"

[[StaticArraysCore]]
git-tree-sha1 = "ec2bd695e905a3c755b33026954b119ea17f2d22"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.3.0"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5950925ff997ed6fb3e985dcce8eb1ba42a0bbe7"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.18"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArraysCore", "Tables"]
git-tree-sha1 = "8c6ac65ec9ab781af05b08ff305ddc727c25f680"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.12"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[Suppressor]]
git-tree-sha1 = "c6ed566db2fe3931292865b966d6d140b7ef32a9"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.1"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[ThreadPools]]
deps = ["Printf", "RecipesBase", "Statistics"]
git-tree-sha1 = "50cb5f85d5646bc1422aa0238aa5bfca99ca9ae7"
uuid = "b189fb0b-2eb5-4ed4-bc0c-d34c51242431"
version = "2.1.1"

[[TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "RecipesBase", "Scratch", "Unicode"]
git-tree-sha1 = "d634a3641062c040fc8a7e2a3ea17661cc159688"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.9.0"

[[ToeplitzMatrices]]
deps = ["AbstractFFTs", "FFTW", "LinearAlgebra", "StatsBase"]
git-tree-sha1 = "39731ea29652305b1eeb7ab0d1085533db307191"
uuid = "c751599d-da0a-543b-9d20-d0a503d91d24"
version = "0.6.3"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "8a75929dcd3c38611db2f8d08546decb514fcadf"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.9"

[[URIs]]
git-tree-sha1 = "e59ecc5a41b000fa94423a578d29290c7266fc10"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─99e0b1b9-d032-41d6-8a1d-742050d5ddeb
# ╟─0a31eff6-d5b0-47c2-8483-042770a77478
# ╠═eeca0160-d37e-11eb-05f4-bd7469319dad
# ╟─0ac4cb74-c027-4677-bd45-672f19035ce4
# ╟─fba8e77d-9141-4f57-afdd-141bb1311a22
# ╟─25b6e79e-37dc-4a67-8b56-fbaf9da41713
# ╟─a4f24d72-fb51-47cc-827d-c9f3cb367998
# ╟─91fa7272-3e83-45c0-b787-87904d4381fb
# ╟─29f3f1a7-c6c1-427a-8781-86ee7dbcd219
# ╟─7e662cb3-8500-4381-ad13-0469f0dbf7eb
# ╟─2cd62a7a-7573-4511-8c81-488de74ad8f5
# ╟─3eaa46e2-5847-4b4d-b48d-fe1142545a8b
# ╟─87a1c372-1e7f-4d00-b897-7739ffdf19e8
# ╟─9771216f-9c98-42c7-98a9-ebc264a5ec0c
# ╟─db94347c-df46-43ee-a7c0-0259c11c859f
# ╟─ef77b7b9-c2ee-4222-a926-5672678aedeb
# ╟─770e4972-6927-4370-8258-25fc101461af
# ╟─9ffc49e0-056b-4029-87ce-2fee216277a5
# ╟─31bc8646-bf07-4a22-9030-196af852251c
# ╟─8313e1ee-7c01-4eda-8d4b-888f04de5c90
# ╠═b6f096f4-9857-4d66-a62d-8dd818c4d9a2
# ╟─7fc20a57-f784-4336-a71f-bb2b0bbdc777
# ╟─55cfd73a-7015-4936-b229-ddcc488d48de
# ╟─59dffcb5-b978-495e-a25e-22a57ae84882
# ╟─ef709d7b-107b-43bb-9026-c7af1efacfdc
# ╟─2e0a22df-6eab-465b-b23b-a9b14a093290
# ╟─18830fc6-d17d-45f9-8e21-cd42ef74ea4e
# ╟─36de0424-f52d-40ca-9c42-dfa05277effa
# ╟─e5910bba-3703-4eb0-9aae-55554193de20
# ╟─095e025d-0344-4e2c-97d9-fa1b230b2469
# ╟─75ba3ecf-6dff-482d-99f3-9650b6597540
# ╟─be3e84ef-842f-43ed-833b-7faee299c904
# ╟─41cee516-36a5-4af3-b157-46fcd1ff6aa5
# ╟─93a84049-7e73-41c4-afd8-0359656c6399
# ╟─95125a71-1691-4009-acb4-5db30cd6b69f
# ╟─507628c6-2c13-4ff5-a73f-c6f30c8e5a8c
# ╟─b1f875af-4046-445e-a2fb-afca07f992ef
# ╟─8918fde2-07d7-4f76-8c5b-03fa04650058
# ╟─c4e41e78-e554-47b3-822f-3516f615bc2c
# ╟─1d8391b6-c28b-4707-9a27-b6152b14be49
# ╠═606948ed-57bd-4f25-8b8a-99133c405bc8
# ╠═c232b1db-ab1c-4c29-b0ee-8fc719b9d2bc
# ╟─79b677e8-98a7-446e-babb-c9f8ccdec52b
# ╟─e45843c7-a0a9-40fd-abff-d150fc1efb67
# ╟─4e21ef20-c5ae-46c7-a512-efb47c48a244
# ╟─46bf582a-ba9a-491f-96f9-fa1ec38698a9
# ╟─53bee8bb-5469-4f99-86f1-a6e9ce0f44f7
# ╟─75885058-267f-4e7e-9cb2-bd5e1f2c87c2
# ╠═456bcf5c-25a2-453d-80ac-3b103a671d93
# ╠═3606489f-e45d-4f06-ab72-2c31def7ed10
# ╠═93d60563-a66e-4a8f-a11b-a81564cd3941
# ╠═b089c56d-9f9e-48c8-bf58-611b58bdcd0c
# ╟─2ca3cba6-90db-4d4b-a225-7ed25093fedd
# ╟─d8b7067c-a4be-4d87-bd79-adb721c61935
# ╟─3bd260fc-de99-4cf9-9453-828df9178b64
# ╠═ec27fb6b-1cbb-46c8-9da6-bf128dd5feae
# ╠═ebfa42bb-0c34-40f1-82ab-dd40d8954eb9
# ╠═a2ef0e9f-8dc9-4892-acbb-3b7b13681c93
# ╠═b1f8acd8-3de8-4512-91e1-47b45780c0ea
# ╟─38a7ad81-944a-438a-ba07-75c24daffdd1
# ╟─113e9d8d-17cc-413e-95ca-0c53dda1f48b
# ╟─ef9fb41f-b815-4e77-8913-82ee14fc3d8c
# ╟─1d3fa479-5844-49ae-bb3c-5a12954cdf60
# ╟─37f7640f-13ab-4d39-a9f1-7a258fe62bf6
# ╟─cb81c6b1-a8ec-473d-bc27-1bdbb46b5a35
# ╠═97b9730f-fffe-4e8d-86e7-553cdd8dcf5b
# ╟─f371ee9e-c72e-4dc0-90d2-b99345edd4d0
# ╟─e46bc8e3-10d7-4637-8d68-9b9ad5f9bf41
# ╠═da986aee-d913-4c23-b214-dd4ae0ede309
# ╠═d27de7be-46b7-4801-b3c3-920ae3c00efe
# ╠═9b1b2070-bda5-43d9-be31-b83829fcd98d
# ╟─59214606-5ad7-4905-b029-43f4271fe5b5
# ╟─a96fe789-13d8-42f3-9c04-07e56637949c
# ╟─f7dba6aa-e1fd-451c-9ab7-e76a483f6959
# ╠═2fa769ef-29d7-440a-89fb-c24980e609ca
# ╠═f0786957-9371-4569-9631-e173b068530e
# ╟─a6837eae-68ac-420b-b47e-a05a32e355a5
# ╠═70edc423-3403-4652-b7a6-e90720543918
# ╠═11b12320-f60e-4a6f-a2fb-b85a573abdd2
# ╟─7387b11f-873b-4af1-b909-bf3051d16e6e
# ╟─5e55cf7d-a076-45e1-8c54-c561be82857d
# ╟─e5ab6234-6221-4ef2-b2b1-05f7e8751928
# ╟─d74faeb0-2338-4f9d-9cae-a240afce8c83
# ╟─fc1ca90d-acbe-4a59-984f-4f9afdc18c42
# ╟─01979f83-69f4-470a-9f47-8d89663b9b9e
# ╟─62a9d6b6-1950-4866-9e2f-86fde888060e
# ╟─6291861f-e52a-4c74-a0bf-fed791d54ef1
# ╟─65a2c53c-3a15-4f48-8cb7-7eb3b3c1b6fb
# ╟─3cf71a1d-0566-4f10-9e71-12f4b9b92eb8
# ╠═9e76aba3-ba12-4212-ad51-2e0b0dd0dee7
# ╟─6c4c01c4-78a0-4a7b-9268-743f8fcf763b
# ╟─02edea51-88dd-4ead-af57-9cf844e70d82
# ╟─19962772-9ddf-4c16-bb2a-bee7abd041e1
# ╟─41a2e257-699a-4c1f-b6de-3441aa0a3cfc
# ╟─0a9c73e1-68c1-4c17-bd48-98673b44d332
# ╠═5dad0cee-8239-43df-b1df-bf47b2d7ec8b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
