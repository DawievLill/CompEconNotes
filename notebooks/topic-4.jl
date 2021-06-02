### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

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
			Pkg.PackageSpec(name="Plots", version="1"), 
			Pkg.PackageSpec(name="LaTeXStrings"), 
			Pkg.PackageSpec(name="Optim"),
			Pkg.PackageSpec(name="OptimTestProblems"),
			Pkg.PackageSpec(name="LineSearches"),
			Pkg.PackageSpec(name="Distributions"), 
			Pkg.PackageSpec(name="IntervalRootFinding"),
			Pkg.PackageSpec(name="Roots"), 
			Pkg.PackageSpec(name="ForwardDiff"),
			Pkg.PackageSpec(name="Convex"),
			Pkg.PackageSpec(name="SCS"),
			Pkg.PackageSpec(name="JuMP"), 
			Pkg.PackageSpec(name="Statistics")
			])
	using Random
	using BenchmarkTools
	using Images
	using PlutoUI
	using HypertextLiteral
	using LinearAlgebra
	using Plots
	using LaTeXStrings
	using Optim
	using OptimTestProblems
	using Distributions
	using IntervalRootFinding
	using Roots
	using ForwardDiff
	using LineSearches
	using Convex
	using SCS
	using JuMP
	using Statistics
end

# ╔═╡ b214155c-17ca-4479-886e-14a09bc1e14c
html"""
<style>
  main {
    max-width: 900px;
  }
</style>
"""

# ╔═╡ 1f4407d0-9b73-11eb-0e91-cd0de83535aa
md" # Optimisation "

# ╔═╡ 16478345-300b-460b-8198-faccaf7740e9
md" We will probably need about two to three sessions to cover the basic optimisation ideas that are relevant for economics. Optimisation is an incredibly important topic in economics, and many related disciplines (such as statistics and machine learning), so it is worthwhile to spend some time exploring the general themes.

Julia has a well-developed ecosystem for optimisation routines and can comfortably used to meet all your optimisation needs. We will speak to the packages relevant for optimisation later in this session. First, we will cover some basics on optimisation and then move to a discussion of several algorithms relevant for optimisation. 

Unfortunately we will only have time to cover a very small subset of the techniques that are relevant to our field. I strongly recommend that you explore this topic in more detail. Some resources will be provided toward the end of the session. The best book on the topic that utilises Julia code is Algorithms for Optimisation, which can be found [here](https://mykel.kochenderfer.com/textbooks/). "

# ╔═╡ 24d0123f-acc4-4656-b564-d7677bd10cf8
md" ## Basics "

# ╔═╡ 9cf3916d-9436-443e-9ed8-75ce7af69d85
md" For this next section I borrow heavily from the excellent notes of **Florian Oswald**, which can be found [here](https://floswald.github.io/NumericalMethods/lecture3/). Let us start with a gentle reminder of optimisation theory, so that we understand what our algorithms are trying to achieve. We will obviously only scratch the surface of optimisation theory, this is a subject in its own right."

# ╔═╡ 4d0b0840-4545-490f-b1f4-84234ef914c1
md" The generic optimisation (minimisation) problem is the following:

$\min _{x \in \mathbb{R}^{n}} f(x) \quad \text{s.t.} \quad x \in \mathscr{X}$

where $x$ is our choice variable, $\mathscr{X}$ is the feasible set, $f$ is the objective function and $x^{*}$ is the solution vector to this problem if it is feasible and minimises $f$."

# ╔═╡ e6b0edda-43ec-4b87-86d8-1ad183d8ab54
md" In economics we often work with constrained optimisation problems where agents need to optimise with regard to some constraint. Firms minimise costs and maximise profits, consumers maximise utility and the social planner maximises social welfare. One such problem relating to utility maximisation given a budget constraint, could be represented as follows:

$\max _{x_{1}, x_{2}} u\left(x_{1}, x_{2}\right) \quad \text{s.t.} \quad p_{1} x_{1}+p_{2} x_{2} \leq y$ 

The constraints define the feasible set. Consider the explicitly defined utility function: 

$\min _{x_{1}, x_{2}}-\exp \left(-\left(x_{1} x_{2}-3 / 2\right)^{2}-\left(x_{2}-3 / 2\right)^{2}\right) \quad \text{s.t.} \quad x_{2} \leq \sqrt{x_{1}}$"

# ╔═╡ 68899349-671d-4166-81c4-4f6fa50f8c46
md" The graph for this constrained optimisation problem is showcased below. The code for the graph is shown. You can view it to gain some insight into plotting in Julia."

# ╔═╡ 9475434f-37bd-413c-85dd-f24f1dd8e504
begin
	gr()
	xx = 0:0.01:3.5 # Range of x values from 0 to 3.5. Incremented by 0.01. 
	f0(x1, x2) = -exp.(-(x1.*x2 - 3/2).^2 - (x2-3/2).^2) # Specification of the function as defined above
	c(z) = sqrt(z) # Constraint function 
	
	p1 = surface(xx, xx, f0, xlab = L"x_1", ylab = L"x_2") # Surface plot
	p2 = contour(xx, xx, f0, lw = 1.5, levels=[collect(0:-0.1:-0.85)...,-0.887,-0.95,-1], lab = L"x_1", ylab = L"x_2") # Contour plot
	
	plot!(p2, c, 0.01, 3.5, label="", lw=2, color=:black, fill=(0,0.5,:blue)) # Constraint
	scatter!(p2, [1.358], [1.165], markersize=5, markercolor=:red, label="Constr. Optimum") # Optimisation point
	plot(p1, p2, size=(900,300)) # Control size of plots
end

# ╔═╡ 90157f5a-f352-4bf7-994f-ee6b1b9ef710
md" #### Optimisation methods "

# ╔═╡ dd6d6b21-5219-4a8e-952b-5446f94d20c1
md" All numerical optimisation methods search through the space of feasible choices, generating guesses that should converge to the true solution. However, methods differ in the information that they use. We have derivative free, comparison methods. Then there are another class of methods that utilise gradients. Thirdly, there are methods that use both the gradient and curvarture about objectives and constraints. 

In order to calculate gradients and curvature we can use numerical, symbolic or automatic differentiation. See our previous lecture on these topics. "

# ╔═╡ c1e1db15-eddb-4535-8d00-c1180a91fdae
md" ### Optimisation packages in Julia "

# ╔═╡ 3b20c100-ff95-4771-957f-634f43d2ccd6
md" Before we move to a discussion on the different algorithms, it might be worthwhile to talk about a few pacakges in Julia that are intended for optimisation. We will mainly focus on the Optim.jl package. "

# ╔═╡ bcf4afe0-e606-4a36-b452-02a7e4f9c7a0
md" #### Optim.jl "

# ╔═╡ be7d0137-a4b9-4625-a629-1b9a69054eea
md" A good package for unconstrained optimisation problems with respect to univariate or multivariate functions is [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl). The default with this package is minimisation, so take that into account when you are considering your optimisation problem.  "

# ╔═╡ 4dbde0ae-88d8-4ed8-baf1-f35b92e6a0e2
md" Let us consider the Rosenbrock function (which we will refer to often in these sessions). The default solver is the Nelder-Mead method, which we will cover in another session. "

# ╔═╡ 442bf711-122b-476f-9695-74e94b8f4df7
begin

	plotly()
	function rosenbrock(t::Vector)
	  return (1.0 - t[1])^2 + 100.0 * (t[2] - t[1]^2)^2
	end
	
	default(size=(600,400), fc=:heat, palette = :Dark2_5)
	t, s = -1.5:0.1:1.5, -1.5:0.1:1.5
	v = Surface((t,s)->rosenbrock([t,s]), t, s)
	surface(t,s,v, linealpha = 0.4)
end

# ╔═╡ ae9a3167-d239-43a5-bedd-f025f3666f47
optimize(rosenbrock, zeros(2)) # The zeros are the initial point for this algorithm to start with. 

# ╔═╡ 67a0a485-bef8-43b0-b8c2-91a0f24fea35
md" We can also use a quasi-Newton method such as L-BFGS. This requires calculation of a gradient, which is done through central finite differencing in the Optim package."

# ╔═╡ 0c2d9cc5-8925-40bd-b101-47a69f1ad2c7
optimize(rosenbrock, zeros(2), LBFGS())

# ╔═╡ 48939170-b144-4035-a8a8-64e2b6caeb64
md" We can specify automatic differentiation as alternative to the default central difference method." 

# ╔═╡ acf3ae5e-9d0a-4fc0-94ec-bac5621be031
optimize(rosenbrock, zeros(2), LBFGS(); autodiff = :forward)

# ╔═╡ 99960cca-a855-4afb-a61a-efefbf53ba90
md" It is worhtwhile to go through the Tutorials in for Optim.jl, they have some great advice about when to use which algorithm. Some other package that are used for optimisation are [JUMP.jl](https://github.com/jump-dev/JuMP.jl), [Roots.jl](https://github.com/JuliaMath/Roots.jl) and [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl/). We will reference these packages as they are required for our purposes. "

# ╔═╡ 38a82b17-13d2-435b-85c8-24c2c619d922
md" # Algorithms for optimisation"

# ╔═╡ 9f559408-6f48-4eca-aa6f-1b3adabd6f06
md" All the algorithms that are used in this section will rely on some form of an iterative procedure. By this we mean that the algorithms will try and improve the value of some objective function over a succession of steps. The way in which the next step is generated is the distinguishing feature of these algorithms. Some algorithms use only objective functions, others include gradients while Hessians are introduced in other methods. "

# ╔═╡ af5cd2a0-eeac-48e3-aff4-a21d2669c75d
md" **Side note**: We will take a look at how to evaluate algorithms in terms of robustness, efficiency and accuracy at a later stage. The best would be if everyone could take a course on algorithms and data structures in the computer science department, but this is perhaps a bit much to ask. "

# ╔═╡ 038a5068-7629-4c91-96f4-081b1e0c6d19
md" We have seen that there are many tools to solve systems of linear equations. However, once the equations become non-linear then things become less obvious. Most of our solution methods for these types of equations work by reducing non-linear equations to sequences of linear equations. We can often approximate our non-linear function by a linear one to get a solution. We can iteratively repeat the process to get better solutions. This is what is known as an iterative algorithm. One example of this is Newton's method, but this requires calculation of derivatives (which can easily be done these days with automatic differentiation).

One of the most basic numerical operations in computational economics is to find the solution to a system of non-linear equations. These non-linear equations can arise in one of two forms. The root finding problem or the non-linear fixed point problem.  The two forms of problems are equivalent, as we specify below. "

# ╔═╡ 0170f06a-79ba-4ba3-bfe0-db1e63bc7b74
md" ## One-dimensional algorithms"

# ╔═╡ 46095078-d68c-4570-b2b5-b8aad3656186
md" The goal of the root finding problem can be stated as follows: 

Given a function $f: \mathbb{R}^{n} \rightarrow \mathbb{R}^n$, find the $x \in \mathbb{R}^{n}$ such that $f(x) = \mathbf{0}$. 

**Example:** If you think about the first model you did in first-year economics, you were solving for the prices and quantities. Clearing markets in a general equilibrium model entails finding the prices such that excess demand functions are zero. 

We also know that fixed point problems are root finding problems, since we can have a non-linear fixed point problem specified as $x = g(x)$, where $g: \mathbb{R}^{n} \rightarrow \mathbb{R}^{n}$ is given and must compute an $n$-vector $x$ called the fixed point of $g$ that satisfies $x = g(x)$. 

This means we can view the root finding problem as a fixed point problem by setting $g(x) = x - f(x)$. Conversely, the fixed point problem can be recast as a root finding problem by letting $f(x) = x - g(x)$.

We encounter root finding problems in **optimisation problems** when the first derivative of the function is a zero. It is also possible to see observe non-linear equations in dynamic optimisation, when we consider collocation methods in the solution of Euler functional equations.  

We start our discussion with one dimensional problems and then move on to higher dimensions. The reason for starting with one dimensional problems is that many multivariate methods reduce to solving sequence of one-dimensional problems. "

# ╔═╡ 6ff4383c-2f61-4e8c-b585-844d0b57543c
md" #### Bracketing methods: Bisection (1D) "

# ╔═╡ c458576d-e3bb-40ea-bb63-0d3690d220d4
md" The simplest version of a root-finding algorithm is the bisection method. This is a gradient free method. The bisection method maintains a bracket $[a, b]$ in which at least one root is known to exist. From the intermediate value theorem we know that we can find a root of the function. 

The **intermediate value theorem** states that if $f$ is a continuous function whose domain contains the interval $[a, b]$, then it takes on any given value between $f(a)$ and $f(b)$ at some point within the interval.

This means that if $f$ is continuous on $[a, b]$, and there is some $y \in [f(a), f(b)]$, then the intermediate value theorem stipulates that there exists at least one $x \in [a, b]$, such that $f(x) = y$. It follows that a bracket $[a, b]$ is guaranteed to contain a zero if $f(a)$ and $f(b)$ have opposite signs. 

This method cuts the bracketed region in half with every iteration. The midpoint is then evaluated, and the new bracket is formed from the midpoint and whichever side that continues to bracket a zero. 
"


# ╔═╡ 51433935-3a0a-42b7-bae4-923006d51bde
md" The general structure of the algorithm is as follows. 

1. Select midpoint of $[a, b] \rightarrow (a + b)/2$
2. Zero must be in lower or upper half. 
3. Check sign of midpoint, if same sign as lower bound a root must be in the right subinterval. 
4. Select midpoint of $[(a + b) / 2, b]$ and continue..."

# ╔═╡ 77b02c2e-b668-4433-9cea-ae7f07b98bea
begin
	function bisection(f, a, b, ε)
		
		if a > b; a, b = b, a; end # ensure a < b
		ya, yb = f(a), f(b)
	
		if ya == 0; b = a; end
		if yb == 0; a = b; end
	
		while b - a > ε
			x = (a + b)/2
			y = f(x)
			if y == 0
				a, b = x, x
			elseif sign(y) == sign(ya)
				a = x
			else
				b = x
			end
		end
		return(a, b)
	end
end

# ╔═╡ 439a4551-13bf-4503-98dd-0f802b3d85dc
md" While we can write our own function, like the one above. It is almost always better to ustilise existing packages. In this case the Roots.jl package will be quite useful. Another package that we can use in the case of multivariate functions is [IntervalRootFinding.jl](https://github.com/JuliaIntervals/IntervalRootFinding.jl/)"

# ╔═╡ bd444293-1a41-4be1-8159-9252dcc6a9db
g(x₁) = exp(x₁) - x₁^4

# ╔═╡ bcd7d97a-8eac-4585-8f38-e712e8038138
plot(g, 8, 9, label = "g")

# ╔═╡ d46b9532-c6b3-4f8d-87dc-03e3a5570656
find_zero(g, (8,9)) # Bisection is the default method.

# ╔═╡ c9a5f2ea-5aec-4214-95b2-ea15d922a1dc
bisection(g, 8, 9, 0.000001) # Result is quite close that of the Roots package. .

# ╔═╡ d3cf3296-7ba5-4a82-a3a8-b8ac83b6e27d
md" Some of the benefits of the bisection method is that it is relatively fast, simple and guaranteed to find a solution if the function satisfies the conditions for the intermediate value theorem. However, this method is difficult to extend to more than one dimension. "

# ╔═╡ 518e561e-9db1-48dd-aa53-f5a55b2c1ca3
md" #### Iterative methods: Newton's method in 1D "

# ╔═╡ 51d7f59a-9958-4401-9ff5-b7874fc48188
md" Newton's method uses information about derivatives by linearising and finds zeros of the newly linearised version. This method is based on the idea of successive linearisation. This method calls for hard non-linear problems to be replaced with a sequence of simpler linear problems, whose solution converges to the solution of the non-linear problem. This problem is normally formulated as a root finding technique, but may be used to solve a fixed point problem by recasting it."

# ╔═╡ e150a318-60ca-428d-9e9c-4c06ef37e231
md" If we have a root approximation $x_k$, we can construct a linear model of $f(x)$ using the formual for the tangent line of a differentiable function (Taylor expansion of the function $f$ around the point $x$): 

$f(x) ≈ f\left(x_{k}\right)+f^{\prime}\left(x_{k}\right)\left(x -x_{k}\right)$

Since $x$ is a zero of the function, we have that,

$0 =  f\left(x_{k}\right)+f^{\prime}\left(x_{k}\right)\left(x -x_{k}\right)$

Our next approximation $f(x_{k+1}) = 0$ means replacing $x = x_{k+1}$, which yields an iteration rule,  

$x_{k+1} \leftarrow x_{k} - \frac{f(x_{k})}{f^{\prime}(x_{k})}$

Starting with an initial estimate $x_1$, this formula produces a sequence of estimates $x_2, x_3, x_4, \ldots$"

# ╔═╡ a103d9a0-332d-4032-9644-de961f36dd8a
md" Let us look at a graphical example to get a feeling for what the algorithm does. "

# ╔═╡ bc9bfd83-b3b2-4f66-aaa3-0f9edf0e651a
straight(x0, y0, x, m) = y0 + m * (x - x0)

# ╔═╡ 866cf80c-26a1-421f-98f1-702bd2de2bb4
function standard_Newton(f, n, x_range, x0, ymin=-10, ymax=10)
    
    f′ = x -> ForwardDiff.derivative(f, x)


	p = plot(f, x_range, lw=3, ylim=(ymin, ymax), legend=:false, size=(400, 300))

	hline!([0.0], c="magenta", lw=3, ls=:dash)
	scatter!([x0], [0], c="green", ann=(x0, -5, L"x_0", 10))

	for i in 1:n

		plot!([x0, x0], [0, f(x0)], c=:gray, alpha=0.5)
		scatter!([x0], [f(x0)], c=:red)
		m = f′(x0)

		plot!(x_range, [straight(x0, f(x0), x, m) for x in x_range], 
			  c=:blue, alpha=0.5, ls=:dash, lw=2)

		x1 = x0 - f(x0) / m

		scatter!([x1], [0], c="green", ann=(x1, -5, L"x_%$i", 10))
		
		x0 = x1

	end

	p |> as_svg


end

# ╔═╡ 0af7eb48-153f-4d57-8691-82403ee3454d
md"""
n = $(@bind n2 Slider(0:10, show_value=true, default=0))
"""

# ╔═╡ 637c4017-a537-4b79-902b-a4bf508e559e
md"""
x₀ = $(@bind x02 Slider(-10:10, show_value=true, default=6))
"""

# ╔═╡ 59372abe-baf4-45f2-b7bd-339ae4dba2bb
let
	gr()
	f(x) = x^2 - 2

	standard_Newton(f, n2, -1:0.01:10, x02, -10, 70)
end

# ╔═╡ f560279a-1c1c-45f9-991d-ceefc575da3d
md" Another example using sliders. "

# ╔═╡ 8d02a306-7ba5-4513-a2e4-103c5637a1ac
md"""
n = $(@bind n Slider(0:10, show_value=true, default=0))
"""

# ╔═╡ 947f1f04-5d95-4047-8f14-1947b3178b30
md"""
x₀ = $(@bind x0 Slider(-10:10, show_value=true, default=6))
"""

# ╔═╡ a6f17a60-2251-427b-b6b3-bacf01927ed0
let
	gr()
	f(x) = 0.2x^3 - 4x + 1
	
	standard_Newton(f, n, -10:0.01:10, x0, -10, 70)
end

# ╔═╡ 3a9f8b1f-bc71-4ae1-93ba-601553d1f4bb
md" One way in which we could code up Newton's method is the following: "

# ╔═╡ 5fa580cd-93fd-4a40-b816-544e7c2acf7e
function newton(f,dfdx,x1)
    # Operating parameters.
    funtol = 100*eps();  xtol = 100*eps();  maxiter = 40;

    x = [x1]
    y = f(x1)
    dx = Inf   # for initial pass below
    k = 1

    while (abs(dx) > xtol) && (abs(y) > funtol) && (k < maxiter)
        dydx = dfdx(x[k])
        dx = -y/dydx            # Newton step
        push!(x,x[k]+dx)        # append new estimate

        k = k+1
        y = f(x[k])
    end

    if k==maxiter
        @warn "Maximum number of iterations reached."
    end

    return x
end

# ╔═╡ 3bcbb714-2b8b-44c7-83c3-568e903a64fe
md" A much easier way is to use automatic differentiation. "

# ╔═╡ 6005e783-2f1b-42a8-9771-cd29ca264e0f
function newton1D(f, x0)
	
	f′(x) = ForwardDiff.derivative(f, x)   
	
	x0 = 37.0  # starting point
	sequence = [x0]
	
	x = x0
	
	for i in 1:10
		x -= f(x) / f′(x)
	end
	
	return x
	
end

# ╔═╡ 57de80f1-b660-4589-b130-6eedd80280f7
newton1D(x -> x^2 - 2, 37.0)

# ╔═╡ 35bab325-8916-49b4-8f46-abd9cd77d4d9
begin
	fx = x -> x^2 - 2;
	dfdx = x -> 2x;
	newton(fx, dfdx, 37.0)
end

# ╔═╡ ceb9c187-d67b-48f2-ade0-46744eda2f4d
sqrt(2)

# ╔═╡ ffd22348-7c3a-4832-a0e4-47bc44dd54be
md" The benfit for Newton's method is that it is even faster than bisection and uses information about the derivative. However, it does not always find a solution and requires the ability to compute a derivative. The secant method is an alternative to Newton's nethod that does not calculate the analytical derivative. We will deal with this quasi-Newton method soon. "

# ╔═╡ 3f62345e-6359-4e2a-89f0-9680798c5f30
md" In terms of convergence of our methods, the bisection method illustrates linear convergence while Newton's method has quadratic convergence. Broadly speaking convergence refers to how fast the solution method converges to the root of the equation. The rate of convergence is the rate of decrease of the bias."

# ╔═╡ a09b11cc-e0f9-4722-b655-b2ebb49b5b83
md" #### Secant method "

# ╔═╡ ff9c9daf-874a-4a81-9250-0d1fa7261a93
md" One of the biggest problems with Newton's method is the fact that you have to calculate $f^{\prime}$. However, this has become easier in recent years with automatic differentiation advances. We can circumvent this problem by making the observation that when a step produces and approximate result we are allowed to carry it out approximately. In our case this means we can use a linear approximation of $f^{\prime}$ in the iteration procedure of Newton's method.

With this method you find the root of the linear approximation through the two most recent root estimates. In other words, 

$f(x) ≈ f\left(x_{k}\right)+\frac{f\left(x_{k}\right)-f\left(x_{k-1}\right)}{x_{k}-x_{k-1}}\left(x_{k}\right)\left(x -x_{k}\right)$

which then means that from our discussion of Newton's method, 

$x_{k+1} \leftarrow x_{k}-\frac{f\left(x_{k}\right)\left(x_{k}-x_{k-1}\right)}{f\left(x_{k}\right)-f\left(x_{k-1}\right)}$

Below is an implementation of the secant method."

# ╔═╡ b21182ad-1806-47b1-9fd4-917308cbfbab
function secant(f,x1,x2)
    # Operating parameters.
    funtol = 100*eps();  xtol = 100*eps();  maxiter = 40;

    x = [x1,x2]
    y1 = f(x1); y2 = 100;
    dx = Inf   # for initial pass below
    k = 2

    while (abs(dx) > xtol) && (abs(y2) > funtol) && (k < maxiter)
        y2 = f(x[k])
        dx = -y2 * (x[k]-x[k-1]) / (y2-y1)   # secant step
        push!(x,x[k]+dx)        # append new estimate

        k = k+1
        y1 = y2    # current f-value becomes the old one next time
    end

    if k==maxiter
        @warn "Maximum number of iterations reached."
    end

    return x
end

# ╔═╡ 32e9ca4d-a60e-4269-b60e-9ac524e1850d
md" If you consider some of our examples above, if the initial value is ill-advised then Newton's method may not converge. One can combine the Newton's method (and its variants) with bracketing to guarantee covergence. The best algorithms use a combination of methods with fast convergence properties and the guarantee on convergence provided by bracketing methods. 

**Brent's method** is an example of such an algorithm and is used in many rootfinding applications. 

Next we consider multidimensional unconstrained optimisation problems. "

# ╔═╡ 93fbb5ea-440f-4095-b3d2-998133930fd6
md" ## Multidimensional unconstrained optimisation"

# ╔═╡ 8a238a48-5a0c-4687-a7e6-125575f74720
md" We move to problems that involve optimisation with mutlivariate functions. We start with local models, which is then followed by first and second order models. As we will see, the general approach for unconstrained optimisation problems is to start at some initial value and work through a series of iterates until we have converged with sufficient accuracy. If our function is smooth we can take advantage of information about the function's shape to figure out where we need to move on the next iteration. "

# ╔═╡ 61d235f5-eefe-43ee-9d8d-59de96cbe539
md" ### Direct methods (gradient free) "

# ╔═╡ 37d802a0-4157-4448-8bc1-ca50d33e16c2
md" The simplest class of derivative free optimisation algorithms are direct methods. They are super simple and often times a good starting point. In essence they simply **compare function values**. These are are also known as heuristic search algorithms. The most popular algorithms include grid search, Nelder-Mead simplex, simulated annealing and particle swarm. We will mention simulated annealing and particle swarm when we talk about methods that instroduce some degree of randomness in the section on stochastic methods.

The basic idea behind these algorithms are then as follows: 

1. Evaluate $f$ at a collection of points.
2. Generate new candidate point $\mathbf{x}^{k}$. Replace the current point with $\mathbf{x}^{k}$ if $f(\mathbf{x}^{k})$ is small enough.
3. Stop when function value stops decreasing. "


# ╔═╡ 41181d84-274d-49c1-99a0-8c3155be8257
md" ##### Random Search "

# ╔═╡ e07068c0-da2c-41b4-8d04-655dd31da954
md" Below is an example of a random search algorithm that implements this general algorithmic idea on the Rosenbrock function "

# ╔═╡ 97d6bbee-7bcb-417b-8626-fe71bd84c773
random_search_url = "https://i.imgur.com/rK1RUKL.gif";

# ╔═╡ 0970a83c-b7f3-4ca3-a29e-eb2e47106026
md""" 
$(Resource(random_search_url, :width => 500))
"""

# ╔═╡ 9a8be602-88b4-4c71-94dc-7a6db978e2f3
md" The code used to generate this random search is provided below."


# ╔═╡ aadfa1a8-22b4-4476-9ffa-4be77bc16d8a
function minrandomsearch(f, x0, npoints; var0=1.0, ftol = 1e-6,
                         vtol = 1e-4, maxiter = 1000,
                         vshrink=0.9, xrange=[-2., 3.],
                         yrange=[-2.,6.], animate=true)
  var = var0     # current variance for search
  oldbest = Inf  # smallest function value
  xbest = x0     # x with smallest function vale
  newbest = f(xbest) 
  iter = 0       # number of iterations
  noimprove = 0  # number of iterations with no improvement
  
  while ((oldbest - newbest > ftol || var > vtol) && iter<=maxiter)
    oldbest = newbest
    x = rand(MvNormal(xbest, var),npoints)

    if animate
      # plot the search points
      p = deepcopy(c)
      scatter!(p, x[1,:], x[2,:], markercolor=:black, markeralpha=0.5, legend=false, xlims=xrange, ylims=yrange)
    end
    
    fval = mapslices(f,x, dims=[1])
    (newbest, i) = findmin(fval)
    if (newbest > oldbest)
      noimprove+=1
      newbest=oldbest
    else
      xbest = x[:,i[2]]
      noimprove = 0
    end
    
    if (noimprove > 10) # shrink var
      var *= vshrink
    end
    
    iter += 1
  end
  if (iter>maxiter)
    info = "Maximum iterations reached"
  else
    info = "Convergence."
  end
  return(newbest, xbest, iter, info)
end

# ╔═╡ e4c81c30-7c6e-4a96-8be8-061093989429
begin
	function banana(a,b)
	  x->(a-x[1])^2+b*(x[2]-x[1]^2)^2
	end
	ff = banana(1.0,1.0)
	
	x_0 = [1.0, 3.0]
	result = minrandomsearch(ff, x_0, 20, var0=0.1, vshrink=0.5, vtol=1e-3 )
end

# ╔═╡ 23e88b42-bfbc-474a-bdc1-3fae1f060829
md" ##### Grid search "

# ╔═╡ 3bab6ddf-8a7c-4774-8e76-cbe51b3f8002
md" With this method we simply compare the objective function at a selection of grid points and pick the highest function value. This method is often a good place to start, but it is quite slow and often requires a large number of grid points to be effective. It is quite robust though and if you evaluate along a large enough number of points you will find a global optimiser! This method is often used in economics since our objective functions may not be nice. Usually a first step in multi-alogorithm approach. General structure of algorithm is as follows

1. Take starting value $x^{(0)}$ and define the region of search, $I = (x^{(0)} - a, x^{(0)} + b)$.
2. Impose a discrete grid consisting of point $x^{(k)}$ where $k = 1, \ldots, n$.
3. Compute $f(x^{(k)})$ for all $k$.
4. Return the maximum of $f(x^{(k)})$ as the result.

Let us consider the function $f(x) = -x^{4} + 2.5x^2 + x + 10$. From the first order conditions we can evaluate the critical points analytically. 

$\begin{aligned} f^{\prime}(x)=-4 x^{3}+5 x+1 &=0 \\-4 x\left(x^{2}-1\right)+x+1 &=0 \\(x+1)\left(-4 x^{2}+4 x+1\right) &=0 \\(x+1)\left(x-\frac{1}{2}-\frac{1}{\sqrt{2}}\right)\left(x-\frac{1}{2}+\frac{1}{\sqrt{2}}\right) &=0 \end{aligned}$"

# ╔═╡ 7c79dfd8-9176-4cd7-bb4b-16468a2f6149
function grid_search(fun, lower_bound = 0, upper_bound =  1, ngrid = 10)
	
	grid = collect(range(lower_bound, upper_bound, length = ngrid))
	func = fun.(grid)
	i = argmax(func)
	return grid[i]
end


# ╔═╡ fb564648-4b03-4303-a392-311731b9c997
f1(x) = -x.^4 .+ (2.5)x.^2 .+ x .+ 10

# ╔═╡ 68bfcba3-517e-4797-8f59-19a7ba47f75a
plot(f1)

# ╔═╡ 3ba55917-2e6d-4eec-906e-b9d664aa9263
critical_values = [-1.0,0.5 - 1/sqrt(2), 0.5 + 1/sqrt(2)]

# ╔═╡ d836b9d1-f3b1-4bee-89ef-b0aeef144279
grid_search(f1, -4, 4, 10)

# ╔═╡ ab60ae00-5b14-4547-9e13-15bdf0e7db56
md" Let us attempt a similar approach with the Rosenbrock function."

# ╔═╡ f372e288-6535-404b-ab18-aa3f475de4dd
# Rosenbrock example from the Optim package

begin
	prob = UnconstrainedProblems.examples["Rosenbrock"]
	ro = prob.f  # Function
	g! = prob.g! # Gradient (analytical)
	h! = prob.h! # Hessian (analytical)
end

# ╔═╡ 2e6c4e41-3262-45f2-a87c-48df3b1b6273
begin
	grid = collect(-1.0:0.1:3);  # grid spacing is important!
	grid2D = [[i;j] for i in grid,j in grid];
	val2D = map(ro,grid2D);
	r = findmin(val2D);
	grid2D[r[2]] # This is the minimiser from the grid search process.  
end

# ╔═╡ 49e664aa-3239-450b-b1fe-3efe5bb650c1
md" Grid search is slow and inaccurate, but it picks out the global optimum every time. A more appropriate example would be:

$f(x)=\left\{\begin{array}{l}\exp (x+3) \text { if } x \in(-\infty,-1] \\ 10 x+13 \text { if } x \in(-1,-0.5] \\ 75 x^{3} \text { if } x \in(-0.5,0.5] \\ 5 \text { if } x \in(0.5,1.5] \\ \log (x-1.5) \text { if } x \in(1.5,+\infty)\end{array}\right.$

Example with many different cases. This function is nasty in that it has non-differentiable kinks, discontinuities, multiple local optima and regions where the function is completely flat. In this case discretisation and grid search might be the only feasible option! "



# ╔═╡ 9d593c53-68a2-48ff-b8c0-611bdc5f24c7
function nasty_func(x)
	
	if x <= -1
		exp(x .+ 3)
	elseif -1 < x <= -0.5
		10x .+ 13
	elseif -0.5 < x <= 0.5
		75(x.^3) 
	elseif 0.5 < x <= 1.5
		5
	else 
		log(x .- 1.5)
	end
end
		

# ╔═╡ d1373376-c76f-40d7-9206-4e205856baa2
plot(-5:0.01:5, nasty_func)

# ╔═╡ 664bca88-60cd-4707-8d12-4e09b1f7b7c2
grid_search(nasty_func, -10, 10, 10) # Use our grid search algorithm from before -- increase grid size for better accuracy.  

# ╔═╡ 59d6e725-b084-4d9f-94ab-f02caed07adf
md" It is often the case that economic models have discontinuities or kinks. Estimation procedures may also require working with piecewise functions. As we have metnioned before, most good algorithms start with a point found by grid search. One could also use the starting value to continue with smaller interval of search, as in adaptive grid search or pattern search in multiple dimensions.  "

# ╔═╡ 5bd99323-8ab3-46f1-8860-15a9cb9876c5
md" ##### Nelder-Mead "

# ╔═╡ 742e6f2f-bdab-41b8-8496-d243db501b01
md" Read the following comment from Lagarias et al. (1999) about the Nelder-Mead algorithm. 

> At present there is no function in any dimension greater than one, for which the original Nelder-Mead algorithm has been proved to converge to a minimizer. Given all the known inefficiencies and failures of the Nelder-Mead algorithm [...], one might wonder why it is used at all, let alone why it is so extraordinarily popular. 

When Nelder-Mead works it is quite fast and some important Matlab functions such as `fmincon` and `fminsearch` use it. It is quite popular, but I am not sure why. Even Judd is very sarcastic about the method in his book. I will not spend any time trying to explain it.

There are many other types of direct methods, but our focus for this session will be more on methods that incorporate information about gradients and Hessians. "

# ╔═╡ 005b33ac-16e0-4813-b22f-53498642d600
md" ### Local descent methods "

# ╔═╡ 9a6e7a48-de02-44dd-8bfb-4ba737f472eb
md" Central question of this section: How can we incrementally improve a design point until some convergence criterion is met?

With local descent methods we are looking for a **local model** that provides some guidance in a certain region of $f$ on where to go next. Local models can be obtained from first- or second-order Taylor approximation. These algorithms are referred to as *descent direction methods*. These methods start with a point $\mathbf{x}^{(1)}$ and generate iterates to converge to a local minimum. 

The steps for local descent are as follows:

1. Check whether $\mathbf{x}^{(k)}$ satisfies termination conditions. If it does terminate, otherwise proceed.  
2. Determine descent direction $\mathbf{d}^{(k)}$ using local information (gradient or Hessian). 
3. Determine step size or learning rate $\alpha^{(k)}$. 
4. Compute next design point according to:

$\mathbf{x}^{(k+1)} \leftarrow \mathbf{x}^{(k)} + \alpha^{(k)}\mathbf{d}^{(k)}$

There are many ways in whcih we can go about determining descent directions and step sizes. Let us consider some strategies for choosing $\alpha$."

# ╔═╡ d31a48dc-05bf-41ce-819a-a57c70b30e65
md" ##### Line search: Finding $\alpha$"

# ╔═╡ ec78165c-1908-4bd8-aabe-44d238902a27
md" One approach for determining the value for $\alpha$ is the line search method. Remember that $\alpha$ tells us how far to move (step length). Here we assume that we have already chosen a specific descent direction $\mathbf{d}$. The line search method helps us determine the value for $\alpha$ using the following minimisation:

$\underset{\alpha}{\operatorname{min}} f(\mathbf{x}+\alpha \mathbf{d})$

We are finding the distance to move, $\alpha$ in direction $\mathbf{d}$ that minimises the objective function $f$. This is a univariate optimisation problem and we can therefore choose from our selection of univariate methods discussed above. Brent-Dekker is most often used to solve this type of problem.

In general it is too costly to calculate the exact value for $\alpha$, so trial values are chosen and the one that generates the lowest value for $f$ is chosen."

# ╔═╡ 8193ff1e-8f99-4426-b0f6-8d08b92aa436
md" If one were to try and code up the line search method this is one approach from the Algorithms for Optimisation textbook."

# ╔═╡ 600476d6-d73c-47e1-a303-6e8a0a72e14b
function bracket_minimum(f, x=0; s=1e-2, k=2.0)
    a, ya = x, f(x)
    b, yb = a + s, f(a + s)
    if yb > ya
        a, b = b, a
        ya, yb = yb, ya
        s = -s
    end
    while true
        c, yc = b + s, f(b + s)
        if yc > yb
            return a < c ? (a, c) : (c, a)
        end
        a, ya, b, yb = b, yb, c, yc
        s *= k
    end
end

# ╔═╡ 46291d46-6af7-4f73-95e1-957232ecc316
function line_search(f, x, d)
    if norm(d) ≈ 0; return x; end; objective = α -> f(x + α*d)
    a, b = bracket_minimum(objective)
    α = minimize(objective, a, b)
    return x + α*d
end

# ╔═╡ aae464a0-9853-43c7-9f12-4b3566ecb30d
md" In the Optim.jl package we find the line search option combined with Newton's method. "

# ╔═╡ d0855f4a-a182-4c12-99ba-598f1ea1d746
# Newton's method with line-search

begin
	algo_hz = Optim.Newton(linesearch = HagerZhang())    # Both Optim.jl and IntervalRootFinding.jl export `Newton`
	res_hz = Optim.optimize(ro, g!, h!, prob.initial_x, method=algo_hz)
end

# ╔═╡ 8a4d3dbe-d241-481b-b946-dd43e2e8fc6f
md" ##### Trust regions: Another way to find $\alpha$. "

# ╔═╡ 32cc3726-241a-4c3f-b9f2-f7d8c8a07c88
md" An alternative to linear search is the trust region method. This method looks at a local model that is believed to be reliable. Limits the step taken by line search and predicts improvement associated with taking a step. If we see that the improvement matches the predicted value then we expand the trust region, otherwise the region contracts. 

With trust region methods we first choose the maximum step size and then the step direction. Once this done, we optimise the step size. The next step is found in this approach by minimising a model of the objective function $\hat{f}$ over a trust region centered at $\mathbf{x}$

Radius of the trust region is expanded and contracted based on how well the model predicts function evaluations. The next design point $\mathbf{x}^{\prime}$ is obtained by solving

$\begin{array}{cc}\underset{\mathbf{x}^{\prime}}{\operatorname{min}} & \hat{f}\left(\mathbf{x}^{\prime}\right) \\ \text { subject to } & \left\|\mathbf{x}-\mathbf{x}^{\prime}\right\| \leq \delta\end{array}$

where the trust region is defined by the positive radius $\delta$ and the vector norm. This is a constrained optimisation problem, something that we will take a look at during another session. "

# ╔═╡ 4805770c-631e-4dd3-8570-857ef7565012
trust_url = "https://i.imgur.com/Mmwb8Nx.jpg";

# ╔═╡ b388bb10-5f68-49fb-ad9a-526987d35244
md""" 
$(Resource(trust_url, :width => 450))
"""

# ╔═╡ f3178cb8-1ac4-486a-8bea-c51d3052f5e3
md" The figure above illustrates the idea behid the trust region method. With this mehtod we constrain the next step to lie within a local region. The trsuted region is expanded or contracted based on the predictive performance of models of the objective function."

# ╔═╡ 9fb15232-7863-4919-a38f-cba715433b44
md" The trust region radius is expanded or contracted based on the local model's predictive performance.These methods compare predicted improvement $\Delta y_{\text {pred }}=f(\mathbf{x})-\hat{f}\left(\mathbf{x}^{\prime}\right)$ to actual improvement $\Delta y_{\mathrm{act}}=f(\mathbf{x})-f\left(\mathbf{x}^{\prime}\right)$:

$\eta=\frac{\text { actual improvement }}{\text { predicted improvement }}=\frac{f(\mathbf{x})-f\left(\mathbf{x}^{\prime}\right)}{f(\mathbf{x})-\hat{f}\left(\mathbf{x}^{\prime}\right)}$

The ratio $\eta$ is close to 1 when the predicted step size matches the actual step size. If the ratio is too small then improvement is considered less than expected and trust region radius is scaled down by factor $\gamma_1 < 1$. The argument runs the other way is the improvement is sufficiently large, and region is scaled up by $\gamma_2 > 1$."

# ╔═╡ 64e694f6-9749-4566-a426-a0ce2f766ff5
md" Below is a graphical representation of the trust region optimisation on the Rosenbrock function. "

# ╔═╡ 5f84ade3-fcd4-43c5-a890-cf7b0b176704
trust_2_url = "https://i.imgur.com/z8y02G2.jpg";

# ╔═╡ f298ae15-eb9f-4506-aca6-9999b71cde49
md""" 
$(Resource(trust_2_url, :width => 400))
"""

# ╔═╡ ab859715-2186-4e5d-8c6c-978140b51f1c
md" If you were to code this up on your own, the algorithm would look something as follows. "

# ╔═╡ 76e3ae4e-de78-484f-8529-51740cf4f659
function solve_trust_region_subproblem(∇f, H, x0, δ)
	x = Variable(length(x0))
	p = Convex.minimize(∇f(x0)⋅(x-x0) + quadform(x-x0, H(x0))/2)
	p.constraints += norm(x-x0) <= δ
	solve!(p, SCSSolver(verbose=false), verbose=false)
	return (x.value, p.optval)
end

# ╔═╡ c01ca185-53f7-401a-a1b7-58b76a19d20f
function trust_region_descent(f, ∇f, H, x, k_max;
	η1=0.25, η2=0.5, γ1=0.5, γ2=2.0, δ=1.0)
	y = f(x)
	for k in 1 : k_max
		x′, y′ = solve_trust_region_subproblem(∇f, H, x, δ)
		r = (y - f(x′)) / (y - y′)
		if r < η1
			δ *= γ1
		else
			x, y = x′, y′
			if r > η2
				δ *= γ2
			end
		end
	end
	return x
end

# ╔═╡ ab954e13-c86e-4d6c-a4f4-9471a1efd15a
md" Trust region methodology is also used in the Optim.jl package together with Newton's method. An example is provided below. "

# ╔═╡ 3bc8fa22-c86a-48b2-8fb7-69f9a680ce73
begin
	# Optim.jl has a TrustRegion for Newton (see below for our discussion on Newton's Method)
	
	NewtonTrustRegion(; initial_delta = 1.0, # The starting trust region radius
	                    delta_hat = 100.0, # The largest allowable trust region radius
	                    eta = 0.1, # When rho is at least eta, accept the step.
	                    rho_lower = 0.25, # When rho is less than rho_lower, shrink the trust region.
	                    rho_upper = 0.75) # When rho is greater than rho_upper, grow the trust region (though no greater than delta_hat).
	
	res = Optim.optimize(ro, g!, h!, prob.initial_x, method=NewtonTrustRegion())
end

# ╔═╡ 07bca67f-d9e0-4bb2-98a2-59a44af7708a
md" In the following sections we need to think about how we will find the descent direction $\mathbf{d}$ for these type of models. "

# ╔═╡ 769e3ed7-1a1c-40da-b5f9-7b06cc7d9ef4
md" #### First-order methods "

# ╔═╡ c07acf72-b729-4364-bfcc-9662ff36489b
md" Here we discuss first order methods to select the appropriate descent direction. "

# ╔═╡ 9676b7d5-5e54-4660-8aac-7c79940478f7
md" ##### Gradient descent "

# ╔═╡ c8482288-04fb-4150-809b-db8d09af0e24
md" An obvious choice for the direction is steepest descent. Steepest descent is the opposite of the gradient. Gradient is defined by:

$\mathbf{g}^{(k)}=\nabla f\left(\mathbf{x}^{(k)}\right)$

and the descent becomes,

$\mathbf{d}^{(k)}=-\frac{\mathbf{g}^{(k)}}{\left\|\mathbf{g}^{(k)}\right\|}$

The method searches along the steepest descent direction at every iteration. Minimising with respect to step size results in a jagged path. In this method each direction is orthogonal to the previous direction. This can easily be proved. 

The benefit of this algorithm is that is only requires the gradient of the function, so no Hessian. However, it can be very slow to converge. There are many first order alternatives to steepest descent that converge much faster. One such an example is the conjugate gradient method, which we will discuss next."

# ╔═╡ 5367cec2-9f7e-48e0-a955-5fd00a69ee87
md" It will later become clear that gradient descent is a type of quasi-Newton solver that takes steps according to,

$x^{(k+1)}=x^{(k)}-P^{-1} \nabla f\left(x^{(k)}\right)$

where $P$ is a positive definite matrix. It is clear to see that the direction for this class of problem $\mathbf{d}^{(k)}$ is $P^{-1} \nabla f\left(x^{(k)}\right)$. If $P$ is the Hessian then we get Newton's method. In gradient descent $P$ is simply an appropriately dimensioned identity matrix, such that we go in the opposite direction of the gradient. No curvature information is used. However, as we stated before, this can be quite slow going if the problem is ill-conditioned. One can use preconditioners to resolve some of these issues. One can also introduce the idea of line search with this method, as with any quasi-Newton method. This means the introduction of a scalar $\alpha$ such that the problem becomes, 

$x^{(k+1)}=x^{(k)}- \alpha P^{-1} \nabla f\left(x^{(k)}\right)$

The value for $\alpha$ needs to be determined by the line search algorithm."

# ╔═╡ 431ddad6-1b19-4d6e-8a43-755673c77c33
GradientDescent(; alphaguess = LineSearches.InitialPrevious(),
                  linesearch = LineSearches.HagerZhang(),
                  P = nothing,
                  precondprep = (P, x) -> nothing);

# ╔═╡ e6077548-5eaa-46ae-abcb-aee09252e01e
md"""
steps = $(@bind ixx1 Slider(1:1000, show_value=true, default=0))
"""

# ╔═╡ 4ab6f872-e504-4753-8a75-72ffb91c3d6c
begin
	plotly()
	xx1 = [-1,1.]
	res2 = optimize(ro, xx1, GradientDescent(), Optim.Options(store_trace=true, extended_trace=true))
	contour(-2.5:0.01:2, -1.5:0.01:2, (xx1,y)->sqrt(ro([xx1, y])), fill=true, color=:deep, legend=false)
	xxtracemat1 = hcat(Optim.x_trace(res2)...)
	plot!(xxtracemat1[1, 1:ixx1], xxtracemat1[2, 1:ixx1], mc = :white, lab="")
	scatter!(xxtracemat1[1:1,ixx1], xxtracemat1[2:2,ixx1], mc=:black, msc=:red, lab="")
	scatter!([1.], [1.], mc=:red, msc=:red, lab="")
	scatter!(xx1[1:1], xx1[2:2], mc=:yellow, msc=:black, label="start", legend=true)
end

# ╔═╡ 6673f432-6368-4455-bf46-dad3cc813901
md" ##### Conjugate gradient "

# ╔═╡ 6053a895-b063-42bb-b30b-9e11597059eb
md" As we have seen in the example above, gradient descent can perform poorly in narrow valeys. The conjugate gradient method tries to overcome this jagged path problem. Combines several methods. Essentially minimises a quadratic form of $f$. Next descent direction uses gradient plus some additional information from the previous step. We won't go into detail on the method as the mathematics becomes a bit complex and the intuition is not obvious. We will simply depict the movement in the value in comparison to gradient descent below. "

# ╔═╡ 515efb25-c1d2-4b52-b2d1-64a1c91a6be7
md"""
steps = $(@bind ix1 Slider(1:40, show_value=true, default=0))
"""

# ╔═╡ 53386bd3-9e26-47c0-86ba-9d5683082523
begin
		plotly()
		xx2 = [-1,1.]
		res3 = optimize(ro, xx2, ConjugateGradient(), Optim.Options(store_trace=true, extended_trace=true))
	    contour(-2.5:0.01:2, -1.5:0.01:2, (xx2,y)->sqrt(ro([xx2, y])), fill=true, color=:deep, legend=false)
	
		xtracemat3 = hcat(Optim.x_trace(res2)...)
		plot!(xtracemat3[1, 1:ix1], xtracemat3[2, 1:ix1], mc = :white,  label="Gradient Descent", legend=true)
		scatter!(xtracemat3[1:1,ix1], xtracemat3[2:2,ix1], mc=:black, msc=:red, label="")	
	
	    xtracemat1 = hcat(Optim.x_trace(res3)...)
	    plot!(xtracemat1[1, 1:ix1], xtracemat1[2, 1:ix1], mc = :white, lab="Conjugate Gradient")
	    scatter!(xtracemat1[1:1,ix1], xtracemat1[2:2,ix1], mc=:black, msc=:red, lab="")
	    scatter!([1.], [1.], mc=:red, msc=:red, lab="")
	    scatter!(xx2[1:1], xx2[2:2], mc=:black, msc=:black, label="start", legend=true, alpha = 0.3)
end

# ╔═╡ 9852d394-32b8-4936-9f73-c5921bb1f31d
md" ##### (L-)BFGS "

# ╔═╡ 509776fd-a055-4f90-9dba-e68f42f61347
md" Technically this method is a first order method, since it only requires gradient information. This quasi-Newton method utilises the same iteration scheme as most quasi-Newton solvers, 

$x^{(k+1)}=x^{(k)}- P^{-1} \nabla f\left(x^{(k)}\right)$

In this case an approximation to the Hessian is built using differences in gradient across iterations. There are two versions of BFGS in the Optim.jl package, namely: BFGS and L-BFGS. The latter method does not use the complete history of the iterative procedure to construct $P$, but only the latest $m$ steps. This method is more suitable for large scale problems as the memory requirement to store relevant vectors grow quickly in large problems. "

# ╔═╡ fd5c5545-74e3-43a5-b0f3-6907ab2efa5f
md" #### Second-order methods "

# ╔═╡ b6a7fe10-8fe3-4a4a-a8c1-1c0e918aeed0
md" In this section we discuss the gold standard with respect to multivariate unconstrained optimisation, Newton's method. Newton's method and its many variations are often the most efficient algorithms for optimisation. These methods update by minimising a second order approximation to the function."

# ╔═╡ 799314f0-7f66-45d5-8798-5227155ce0bc
md" ##### Newton's method "

# ╔═╡ a0ce61ac-2478-46ce-b430-9ad5a67888d2
md" Main advantage of this method is that it has quadratic convergence near the local optimum. However, one needs to provide the Hessian with this method, which could be quite costly to do. This method for optimisation consists of applying Newton's method for solving systems of equations, where the equations are now first order conditions. This is a root-finding problem on the first order conditions. This equivalent to saying that the gradient should equal the zero vector.

In the univariate case quadratic optimisation about the point $x^{(k)}$ came from the second-order Taylor expansion 

$q(x)=f\left(x^{(k)}\right)+\left(x-x^{(k)}\right) f^{\prime}\left(x^{(k)}\right)+\frac{\left(x-x^{(k)}\right)^{2}}{2} f^{\prime \prime}\left(x^{(k)}\right)$

Setting the derivative to zero and solving for the root gives us the update equation for Newton's method:

$\begin{aligned} \frac{\partial}{\partial x} q(x) &=f^{\prime}\left(x^{(k)}\right)+\left(x-x^{(k)}\right) f^{\prime \prime}\left(x^{(k)}\right)=0 \\ x^{(k+1)} &=x^{(k)}-\frac{f^{\prime}\left(x^{(k)}\right)}{f^{\prime \prime}\left(x^{(k)}\right)} \end{aligned}$

This can be extended to multivariate optimisation, with Taylor expansion at $\mathbf{x}^{(k)}$ being

$f(\mathbf{x}) \approx q(\mathbf{x})=f\left(\mathbf{x}^{(k)}\right)+\left(\mathbf{g}^{(k)}\right)^{\top}\left(\mathbf{x}-\mathbf{x}^{(k)}\right)+\frac{1}{2}\left(\mathbf{x}-\mathbf{x}^{(k)}\right)^{\top} \mathbf{H}^{(k)}\left(\mathbf{x}-\mathbf{x}^{(k)}\right)$

Evaluate the gradient and set it to zero:

$\nabla q\left(\mathbf{x}^{(k)}\right)=\mathbf{g}^{(k)}+\mathbf{H}^{(k)}\left(\mathbf{x}-\mathbf{x}^{(k)}\right)=\mathbf{0}$

Solve for the next iterate and obtain Newton's method in multivariate form, 

$\mathbf{x}^{(k+1)}=\mathbf{x}^{(k)}-\left(\mathbf{H}^{(k)}\right)^{-1} \mathbf{g}^{(k)}$

We can also write this more generally as, 

$\mathbf{x}^{(k+1)}=\mathbf{x}^{(k)}- P^{-1} \nabla f\left(\mathbf{x}^{(k)}\right)$

where $P$ is the Hessian and $\nabla f\left(x^{(k)}\right)$ is the gradient. Newton's method can also be used to supply a descent direction to line search. The descent direction is,

$\mathbf{d}^{(k)}=-\left(\mathbf{H}^{(k)}\right)^{-1} \mathbf{g}^{(k)}$

Below we provide an example of the different speeds of convergence of gradient descent, conjugate gradient and Newton's method. 
"

# ╔═╡ 978de32a-7975-4d85-8c3e-36e7c2efdd83
md"""
steps = $(@bind ix2 Slider(1:22, show_value=true, default=0))
"""

# ╔═╡ 4d1e8b43-1d74-4100-9eaf-c3923c0e269e
begin
	plotly()
	xx3 = [-1., 1.]
	res4 = optimize(ro, g!, h!, xx3, Optim.Newton(), Optim.Options(store_trace=true, extended_trace=true))
	contour(-2.5:0.01:2, -1.5:0.01:2, (xx3,y)->sqrt(ro([xx3, y])), fill=true, color=:deep, legend=false)
	
	xxtracemat2 = hcat(Optim.x_trace(res2)...)
	plot!(xxtracemat2[1, 1:ix2], xxtracemat2[2, 1:ix2], mc = :white,  label="Gradient Descent", legend=true)
	scatter!(xxtracemat2[1:1,ix2], xxtracemat2[2:2,ix2], mc=:black, msc=:red, label="")
	
	xxtracemat3 = hcat(Optim.x_trace(res3)...)
    plot!(xxtracemat3[1, 1:ix2], xxtracemat3[2, 1:ix2], mc = :white,  label="Conjugate Gradient", legend=true)
    scatter!(xxtracemat3[1:1,ix2], xxtracemat3[2:2,ix2], mc=:red, msc=:black, label="")
	
	xtracemat2 = hcat(Optim.x_trace(res4)...)
	plot!(xtracemat2[1, 1:ix2], xtracemat2[2, 1:ix2], c=:blue, label="Newton")
	scatter!(xtracemat2[1:1,ix2], xtracemat2[2:2,ix2], mc=:black, msc=:blue, label="")
	scatter!([1.], [1.], mc=:red, msc=:red, label="")
	scatter!(xx3[1:1], xx3[2:2], mc=:black, msc=:black, label="start", alpha = 0.3)
	
end

# ╔═╡ 0a37a8a2-b56f-4c30-9abd-5b4472307388
md" Newton's method tends to take relatively few iterations to converge in the case of well-behaved functions. The biggest concern is calculation of Hessians and their inverses. Another area of concern is where functions are not well approximated by their second expansions. Normally we deal with this problem by combining Newton's method with trust regions or line searches."

# ╔═╡ 41babb68-bf77-4679-a3db-dd85c207ee68
@benchmark optimize(ro, [0.0, 0.0], Optim.Newton(),Optim.Options(show_trace=false))

# ╔═╡ 35ec71a1-eb61-4e31-a3bb-6e48d790496f
@benchmark optimize(ro, g!, h!,  [-1.0, 3.0], BFGS())

# ╔═╡ e1696385-ba74-4c5e-85c6-2bf410dc3c70
@benchmark optimize(ro, g!, h!,  [0.0, 0.0], LBFGS())

# ╔═╡ 9715ec3c-0477-4cc0-8a58-b76e8f886895
md" #### Stochastic methods "

# ╔═╡ 5c141a51-2cbf-4312-85df-34042039eba1
md" Here we discuss some method that incorporate dimensions of randomness. We will focus on stochastic gradient descent, simulated annealing and particle swarm methods. "

# ╔═╡ 1133f64e-b997-42b1-8603-caca664baadb
md" ##### Stochastic gradient descent "

# ╔═╡ 7b05748a-6f28-4360-af3b-50c75e13ec64
md" We know that methods like gradient descent can get stuck at local minima. One way to resolve this issue is to apply a random shock to the value of the descent direction. One could modify the gradient from before in the following fashion, 

$\mathbf{x}^{(k+1)} \leftarrow \mathbf{x}^{(k)}+\alpha^{k} \mathbf{g}^{(k)}+\boldsymbol{\varepsilon}^{(k)}$

where $\boldsymbol{\varepsilon}^{(k)} \sim N\left(0, \sigma_{k}^{2}\right)$ is decreasing in $k$. This method is quite popular in training neural networks in deep learning. "

# ╔═╡ 10cd58da-6f5f-4d3b-961f-5c0b7d21ecba
md" ##### Simulated annealing "

# ╔═╡ a501dd90-021c-47e2-92de-9e5795958930
md" In this method one specifies a **temperature** that controls the degree of randomness. If the temperature is high then the search is alowed to jump around across a wide range of values. Start with an initial high temperature in order to escape potential local minima. "

# ╔═╡ def9105a-5eac-49b7-a451-729c94e3baa6
md" ##### Particle swarm "

# ╔═╡ f25f23c0-1436-4744-a05c-a92b10ad25b9
md" ### Multidimensional constrained optimisation"

# ╔═╡ Cell order:
# ╟─f4226cfe-ee06-4c72-9615-fc4aedfd045c
# ╟─b214155c-17ca-4479-886e-14a09bc1e14c
# ╟─1f4407d0-9b73-11eb-0e91-cd0de83535aa
# ╟─16478345-300b-460b-8198-faccaf7740e9
# ╟─24d0123f-acc4-4656-b564-d7677bd10cf8
# ╟─9cf3916d-9436-443e-9ed8-75ce7af69d85
# ╟─4d0b0840-4545-490f-b1f4-84234ef914c1
# ╟─e6b0edda-43ec-4b87-86d8-1ad183d8ab54
# ╟─68899349-671d-4166-81c4-4f6fa50f8c46
# ╠═9475434f-37bd-413c-85dd-f24f1dd8e504
# ╟─90157f5a-f352-4bf7-994f-ee6b1b9ef710
# ╟─dd6d6b21-5219-4a8e-952b-5446f94d20c1
# ╟─c1e1db15-eddb-4535-8d00-c1180a91fdae
# ╟─3b20c100-ff95-4771-957f-634f43d2ccd6
# ╟─bcf4afe0-e606-4a36-b452-02a7e4f9c7a0
# ╟─be7d0137-a4b9-4625-a629-1b9a69054eea
# ╟─4dbde0ae-88d8-4ed8-baf1-f35b92e6a0e2
# ╠═442bf711-122b-476f-9695-74e94b8f4df7
# ╠═ae9a3167-d239-43a5-bedd-f025f3666f47
# ╟─67a0a485-bef8-43b0-b8c2-91a0f24fea35
# ╠═0c2d9cc5-8925-40bd-b101-47a69f1ad2c7
# ╟─48939170-b144-4035-a8a8-64e2b6caeb64
# ╠═acf3ae5e-9d0a-4fc0-94ec-bac5621be031
# ╟─99960cca-a855-4afb-a61a-efefbf53ba90
# ╟─38a82b17-13d2-435b-85c8-24c2c619d922
# ╟─9f559408-6f48-4eca-aa6f-1b3adabd6f06
# ╟─af5cd2a0-eeac-48e3-aff4-a21d2669c75d
# ╟─038a5068-7629-4c91-96f4-081b1e0c6d19
# ╟─0170f06a-79ba-4ba3-bfe0-db1e63bc7b74
# ╟─46095078-d68c-4570-b2b5-b8aad3656186
# ╟─6ff4383c-2f61-4e8c-b585-844d0b57543c
# ╟─c458576d-e3bb-40ea-bb63-0d3690d220d4
# ╟─51433935-3a0a-42b7-bae4-923006d51bde
# ╠═77b02c2e-b668-4433-9cea-ae7f07b98bea
# ╟─439a4551-13bf-4503-98dd-0f802b3d85dc
# ╠═bd444293-1a41-4be1-8159-9252dcc6a9db
# ╠═bcd7d97a-8eac-4585-8f38-e712e8038138
# ╠═d46b9532-c6b3-4f8d-87dc-03e3a5570656
# ╠═c9a5f2ea-5aec-4214-95b2-ea15d922a1dc
# ╟─d3cf3296-7ba5-4a82-a3a8-b8ac83b6e27d
# ╟─518e561e-9db1-48dd-aa53-f5a55b2c1ca3
# ╟─51d7f59a-9958-4401-9ff5-b7874fc48188
# ╟─e150a318-60ca-428d-9e9c-4c06ef37e231
# ╟─a103d9a0-332d-4032-9644-de961f36dd8a
# ╟─bc9bfd83-b3b2-4f66-aaa3-0f9edf0e651a
# ╟─866cf80c-26a1-421f-98f1-702bd2de2bb4
# ╟─0af7eb48-153f-4d57-8691-82403ee3454d
# ╟─637c4017-a537-4b79-902b-a4bf508e559e
# ╟─59372abe-baf4-45f2-b7bd-339ae4dba2bb
# ╟─f560279a-1c1c-45f9-991d-ceefc575da3d
# ╟─8d02a306-7ba5-4513-a2e4-103c5637a1ac
# ╟─947f1f04-5d95-4047-8f14-1947b3178b30
# ╟─a6f17a60-2251-427b-b6b3-bacf01927ed0
# ╟─3a9f8b1f-bc71-4ae1-93ba-601553d1f4bb
# ╠═5fa580cd-93fd-4a40-b816-544e7c2acf7e
# ╟─3bcbb714-2b8b-44c7-83c3-568e903a64fe
# ╠═6005e783-2f1b-42a8-9771-cd29ca264e0f
# ╠═57de80f1-b660-4589-b130-6eedd80280f7
# ╠═35bab325-8916-49b4-8f46-abd9cd77d4d9
# ╠═ceb9c187-d67b-48f2-ade0-46744eda2f4d
# ╟─ffd22348-7c3a-4832-a0e4-47bc44dd54be
# ╟─3f62345e-6359-4e2a-89f0-9680798c5f30
# ╟─a09b11cc-e0f9-4722-b655-b2ebb49b5b83
# ╟─ff9c9daf-874a-4a81-9250-0d1fa7261a93
# ╠═b21182ad-1806-47b1-9fd4-917308cbfbab
# ╟─32e9ca4d-a60e-4269-b60e-9ac524e1850d
# ╟─93fbb5ea-440f-4095-b3d2-998133930fd6
# ╟─8a238a48-5a0c-4687-a7e6-125575f74720
# ╟─61d235f5-eefe-43ee-9d8d-59de96cbe539
# ╟─37d802a0-4157-4448-8bc1-ca50d33e16c2
# ╟─41181d84-274d-49c1-99a0-8c3155be8257
# ╟─e07068c0-da2c-41b4-8d04-655dd31da954
# ╟─97d6bbee-7bcb-417b-8626-fe71bd84c773
# ╟─0970a83c-b7f3-4ca3-a29e-eb2e47106026
# ╟─9a8be602-88b4-4c71-94dc-7a6db978e2f3
# ╠═aadfa1a8-22b4-4476-9ffa-4be77bc16d8a
# ╠═e4c81c30-7c6e-4a96-8be8-061093989429
# ╟─23e88b42-bfbc-474a-bdc1-3fae1f060829
# ╟─3bab6ddf-8a7c-4774-8e76-cbe51b3f8002
# ╠═7c79dfd8-9176-4cd7-bb4b-16468a2f6149
# ╠═fb564648-4b03-4303-a392-311731b9c997
# ╠═68bfcba3-517e-4797-8f59-19a7ba47f75a
# ╠═3ba55917-2e6d-4eec-906e-b9d664aa9263
# ╠═d836b9d1-f3b1-4bee-89ef-b0aeef144279
# ╟─ab60ae00-5b14-4547-9e13-15bdf0e7db56
# ╠═f372e288-6535-404b-ab18-aa3f475de4dd
# ╠═2e6c4e41-3262-45f2-a87c-48df3b1b6273
# ╟─49e664aa-3239-450b-b1fe-3efe5bb650c1
# ╠═9d593c53-68a2-48ff-b8c0-611bdc5f24c7
# ╠═d1373376-c76f-40d7-9206-4e205856baa2
# ╠═664bca88-60cd-4707-8d12-4e09b1f7b7c2
# ╟─59d6e725-b084-4d9f-94ab-f02caed07adf
# ╟─5bd99323-8ab3-46f1-8860-15a9cb9876c5
# ╟─742e6f2f-bdab-41b8-8496-d243db501b01
# ╟─005b33ac-16e0-4813-b22f-53498642d600
# ╟─9a6e7a48-de02-44dd-8bfb-4ba737f472eb
# ╟─d31a48dc-05bf-41ce-819a-a57c70b30e65
# ╟─ec78165c-1908-4bd8-aabe-44d238902a27
# ╟─8193ff1e-8f99-4426-b0f6-8d08b92aa436
# ╠═600476d6-d73c-47e1-a303-6e8a0a72e14b
# ╠═46291d46-6af7-4f73-95e1-957232ecc316
# ╟─aae464a0-9853-43c7-9f12-4b3566ecb30d
# ╠═d0855f4a-a182-4c12-99ba-598f1ea1d746
# ╟─8a4d3dbe-d241-481b-b946-dd43e2e8fc6f
# ╟─32cc3726-241a-4c3f-b9f2-f7d8c8a07c88
# ╟─4805770c-631e-4dd3-8570-857ef7565012
# ╟─b388bb10-5f68-49fb-ad9a-526987d35244
# ╟─f3178cb8-1ac4-486a-8bea-c51d3052f5e3
# ╟─9fb15232-7863-4919-a38f-cba715433b44
# ╟─64e694f6-9749-4566-a426-a0ce2f766ff5
# ╟─5f84ade3-fcd4-43c5-a890-cf7b0b176704
# ╟─f298ae15-eb9f-4506-aca6-9999b71cde49
# ╟─ab859715-2186-4e5d-8c6c-978140b51f1c
# ╠═c01ca185-53f7-401a-a1b7-58b76a19d20f
# ╠═76e3ae4e-de78-484f-8529-51740cf4f659
# ╟─ab954e13-c86e-4d6c-a4f4-9471a1efd15a
# ╠═3bc8fa22-c86a-48b2-8fb7-69f9a680ce73
# ╟─07bca67f-d9e0-4bb2-98a2-59a44af7708a
# ╟─769e3ed7-1a1c-40da-b5f9-7b06cc7d9ef4
# ╟─c07acf72-b729-4364-bfcc-9662ff36489b
# ╟─9676b7d5-5e54-4660-8aac-7c79940478f7
# ╟─c8482288-04fb-4150-809b-db8d09af0e24
# ╟─5367cec2-9f7e-48e0-a955-5fd00a69ee87
# ╟─431ddad6-1b19-4d6e-8a43-755673c77c33
# ╟─e6077548-5eaa-46ae-abcb-aee09252e01e
# ╠═4ab6f872-e504-4753-8a75-72ffb91c3d6c
# ╟─6673f432-6368-4455-bf46-dad3cc813901
# ╟─6053a895-b063-42bb-b30b-9e11597059eb
# ╟─515efb25-c1d2-4b52-b2d1-64a1c91a6be7
# ╠═53386bd3-9e26-47c0-86ba-9d5683082523
# ╟─9852d394-32b8-4936-9f73-c5921bb1f31d
# ╟─509776fd-a055-4f90-9dba-e68f42f61347
# ╟─fd5c5545-74e3-43a5-b0f3-6907ab2efa5f
# ╟─b6a7fe10-8fe3-4a4a-a8c1-1c0e918aeed0
# ╟─799314f0-7f66-45d5-8798-5227155ce0bc
# ╟─a0ce61ac-2478-46ce-b430-9ad5a67888d2
# ╟─978de32a-7975-4d85-8c3e-36e7c2efdd83
# ╠═4d1e8b43-1d74-4100-9eaf-c3923c0e269e
# ╟─0a37a8a2-b56f-4c30-9abd-5b4472307388
# ╠═41babb68-bf77-4679-a3db-dd85c207ee68
# ╠═35ec71a1-eb61-4e31-a3bb-6e48d790496f
# ╠═e1696385-ba74-4c5e-85c6-2bf410dc3c70
# ╟─9715ec3c-0477-4cc0-8a58-b76e8f886895
# ╟─5c141a51-2cbf-4312-85df-34042039eba1
# ╟─1133f64e-b997-42b1-8603-caca664baadb
# ╟─7b05748a-6f28-4360-af3b-50c75e13ec64
# ╟─10cd58da-6f5f-4d3b-961f-5c0b7d21ecba
# ╟─a501dd90-021c-47e2-92de-9e5795958930
# ╟─def9105a-5eac-49b7-a451-729c94e3baa6
# ╟─f25f23c0-1436-4744-a05c-a92b10ad25b9
