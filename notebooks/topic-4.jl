### A Pluto.jl notebook ###
# v0.14.5

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
			Pkg.PackageSpec(name="Plots"), 
			Pkg.PackageSpec(name="LaTeXStrings"), 
			Pkg.PackageSpec(name="Optim"),
			Pkg.PackageSpec(name="Distributions"), 
			Pkg.PackageSpec(name="IntervalRootFinding"),
			Pkg.PackageSpec(name="Roots"), 
			Pkg.PackageSpec(name="ForwardDiff"),
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
	using Distributions
	using IntervalRootFinding
	using Roots
	using ForwardDiff
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
md" For this next section I borrow heavily from the excellent notes of **Florian Oswald**, which can be found [here](https://floswald.github.io/NumericalMethods/lecture3/). Let us start with a gentle reminder of optimisation theory, so that we understand what our algorithms are trying to achieve. We will obviously only scratch the surface of optimisation theory, this is a subject in its own right."

# ╔═╡ 4d0b0840-4545-490f-b1f4-84234ef914c1
md" The generic optimisation (minimisation) problem is the following:

$\min _{x \in \mathbb{R}^{n}} f(x) \quad \text{s.t.} \quad x \in \mathscr{X}$

where $x$ is our choice variable, $\mathscr{X}$ is the feasible set, $f$ is the objective function and $x^{*}$ is the solution vector to this problem if it is feasible and minimises $f$."

# ╔═╡ e6b0edda-43ec-4b87-86d8-1ad183d8ab54
md" In economics we often work with constrained optimisation problems such as the following:

$\max _{x_{1}, x_{2}} u\left(x_{1}, x_{2}\right) \quad \text{s.t.} \quad p_{1} x_{1}+p_{2} x_{2} \leq y$ 

The constraints define the feasible set. Let us consider the following example with explicitly defined utility function: 

$\min _{x_{1}, x_{2}}-\exp \left(-\left(x_{1} x_{2}-3 / 2\right)^{2}-\left(x_{2}-3 / 2\right)^{2}\right) \quad \text{s.t.} \quad x_{2} \leq \sqrt{x_{1}}$"

# ╔═╡ 68899349-671d-4166-81c4-4f6fa50f8c46
md" The graph for this constrained optimisation problem is showcased below. The code for the graph is shown. You can view it to gain some insight into plotting in Julia."

# ╔═╡ 9475434f-37bd-413c-85dd-f24f1dd8e504
begin
	
	gr()
	
	xx = 0:0.01:3.5 # Domain of x values from 0 to 3.5. Incremented by 0.01. 
	f0(x1, x2) = -exp.(-(x1.*x2 - 3/2).^2 - (x2-3/2).^2) # Specification of the function as defined above
	c(z) = sqrt(z) # Constraint function 
	
	p1 = surface(xx, xx, f0, xlab = L"x_1", ylab = L"x_2") # Surface plot
	p2 = contour(xx, xx, f0, lw = 1.5, levels=[collect(0:-0.1:-0.85)...,-0.887,-0.95,-1], lab = L"x_1", ylab = L"x_2") # Contour plot
	
	plot!(p2, c, 0.01, 3.5, label="", lw=2, color=:black, fill=(0,0.5,:blue)) # Constraint
	scatter!(p2, [1.358], [1.165], markersize=5, markercolor=:red, label="Constr. Optimum") # Optimisation point
	plot(p1, p2, size=(900,300)) # Control size of plots
end

# ╔═╡ 90157f5a-f352-4bf7-994f-ee6b1b9ef710
md" #### Critical points "

# ╔═╡ dd6d6b21-5219-4a8e-952b-5446f94d20c1
md" Univariate functions can have several critical points, those where the derivative is zero. We would like to find the global minimum, but this is not always easy to do. Most of the time we will only be finding the **local minimum**.

We have necesarry and sufficient conditions for optimisation. Here we point out the necesarry conditions. In terms of the minimisation problem, for the univariate case we know that the first order derivative needs to be equal to zero, while the second order derivative is required to be non-negative for local optima. For the multivariate case we turn to conditions on gradients instead of derivatives. In this case the graadient is zero and Hessian is positive semidefinite for a local minima.

In order to calculate derivates and gradients we can use numerical, symbolic or automatic differentiation. See our previous lecture on these topics. "

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
	function rosenbrock(t::Vector)
	  return (1.0 - t[1])^2 + 100.0 * (t[2] - t[1]^2)^2
	end
	
	default(size=(600,400), fc=:heat)
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
md" It is worhtwhile to go through the Tutorials in for Optim.jl, they have some great advice about when to use which algorithm. Some other package that are used for optimisation are [JUMP.jl](https://github.com/jump-dev/JuMP.jl), [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl), [Roots.jl](https://github.com/JuliaMath/Roots.jl) and [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl/). We will reference these packages as they are required for our purposes. "

# ╔═╡ 38a82b17-13d2-435b-85c8-24c2c619d922
md" ## Algorithms for optimisation"

# ╔═╡ 9f559408-6f48-4eca-aa6f-1b3adabd6f06
md" All the algorithms that are used in this section will rely on some form of an iterative procedure. By this we mean that the algorithms will try and improve the value of some objective function over a succession of steps. The way in which the next step is generated is the distinguishing feature of these algorithms. Some algorithms use only objective functions, others include gradients while Hessians are introduced in other methods. "

# ╔═╡ af5cd2a0-eeac-48e3-aff4-a21d2669c75d
md" **Side note**: We will take a look at how to evaluate algorithms in terms of robustness, efficiency and accuracy at a later stage. The best would be if everyone could take a course on algorithms and data structures in the computer science department, but this is perhaps a bit much to ask. "

# ╔═╡ 0170f06a-79ba-4ba3-bfe0-db1e63bc7b74
md" #### Root finding algorithms"

# ╔═╡ 46095078-d68c-4570-b2b5-b8aad3656186
md" The goal of the root finding problem can be stated as follows: 

Given a function $f: \mathbb{R}^{n} \rightarrow \mathbb{R}^n$, find the $x \in \mathbb{R}^{n}$ such that $f(x) = \mathbf{0}$. 

**Example:** If you think about the first model you did in first-year economics, you were solving for the prices and quantities. Clearing markets in a general equilibrium model entails finding the prices such that excess demand functions are zero. 

We also know that fixed point problems are root finding problems, so these seem quite important. This means if we are given an equation $g(x) = h(x)$ we can define $f = g -h$ and find the root of $f$

Unlike with linear problems from our earlier discussion, the root will often not be produced in a finite number of steps. So we will need a sequence of approximations that converge to the root, which stops when some member of the sequence is `good enough`. "

# ╔═╡ 6ff4383c-2f61-4e8c-b585-844d0b57543c
md" ##### Bracketing methods: Bisection "

# ╔═╡ c458576d-e3bb-40ea-bb63-0d3690d220d4
md" The simplest version of a root-finding algorithm is the bisection method. This is a gradient free method. The bisection method maintains a bracket $[a, b]$ in which at least one root is known to exist. From the intermediate value theorem we know that we can find a root of the function. 

The intermediate value theorem states that if $f$ is a continuous function whose domain contains the interval $[a, b]$, then it takes on any given value between $f(a)$ and $f(b)$ at some point within the interval.

This means that if $f$ is continuous on $[a, b]$, and there is some $y \in [f(a), f(b)]$, then the intermediate value theorem stipulates that there exists at least one $x \in [a, b]$, such that $f(x) = y$. It follows that a bracket $[a, b]$ is gauranteed to contain a zero if $f(a)$ and $f(b)$ have opposite signs. 

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
bisection(g, 8, 9, 0.1) # Result is quite close that of the Roots package. .

# ╔═╡ d3cf3296-7ba5-4a82-a3a8-b8ac83b6e27d
md" Some of the benefits of the bisection method is that it is relatively fast, simple and guaranteed to find a solution if the function satisfies the conditions for the intermediate value theorem. However, this method is difficult to extend to more than one dimension. "

# ╔═╡ 518e561e-9db1-48dd-aa53-f5a55b2c1ca3
md" ##### Iterative methods: Newton's method in 1D "

# ╔═╡ 51d7f59a-9958-4401-9ff5-b7874fc48188
md" Newton's method uses information about derivatives by linearising and finds zeros of the newly linearised version. Before moving on to the math, let us look at the idea of Newton's method. The following graph provides a function, of which we want to find the zero / root."

# ╔═╡ fd5842db-c6b5-44d5-b60b-eadb51913247
begin
	h = x -> x*exp(x) - 2
	
	plot(h,0,1.5,label="function",grid=:y,xlabel="x",ylabel="y",legend=:topleft)
end

# ╔═╡ 3f1d4727-0f5e-4513-abb9-945eed7f9874
md" We can clearly see that the root is somewhere near $x = 1$, so we make an initial guess at that point." 

# ╔═╡ d2ca4b6c-b82a-4ca7-8b7f-9274ace2c628
begin
	x1 = 1
	h1 = h(x1)
	scatter!([x1],[h1],label="initial point")
end

# ╔═╡ 7e9a56d5-b169-4cfa-afc3-276743615327
md" Next, we compute the tangent line at the point $(x_1, h(x_1))$ using the derivative."

# ╔═╡ 2a1e1787-d850-4ef0-b7c7-3682d0cc8881
begin
	dhdx = x -> exp(x)*(x+1)
	slope1 = dhdx(x1)
	tangent1 = x -> h1 + slope1*(x-x1)
	
	plot!(tangent1,0,1.5,l=:dash,label="tangent line",ylim=[-2,4])
end

# ╔═╡ 62c6879b-e609-45d6-affd-c82e69356919
md" This does not provide the root of $f$, so we settle for finding the root of the tangent line approximation, which is trivial. Call this $x_2$, our next approximation to the root."

# ╔═╡ 91ca8fdc-d3e7-4d9a-8613-f67505805b21
begin
	x2 = x1 - h1/slope1
	scatter!([x2],[0],label="tangent root")
end

# ╔═╡ 812e6f4b-c266-4e29-a5f8-465f46c1b1d2
x2

# ╔═╡ 6a50c48b-e41c-4028-a05f-e4b001dcd786
h2 = h(x2)

# ╔═╡ 749d1d8a-4dd6-4509-81c0-a880d95cb8e1
md" The residual (value of $h$) is smaller than before, but not zero. So we repeat the process with a new tangent line based on the latest point on the curve."

# ╔═╡ ff385d11-7916-4a68-b4d0-2ade41f252e0
begin
	begin
		plot(h,0.8,0.9,label="function",
		    xlabel="x", ylabel="y", title="Second iteration", legend=:topleft)
		
		scatter!([x2],[h2],label="starting point")
		
		slope2 = dhdx(x2)
		tangent2 = x -> h2 + slope2*(x-x2)
		plot!(tangent2,0.8,0.9,l=:dash,label="tangent line")
		
		x3 = x2 - h2/slope2
		scatter!([x3],[0],label="tangent root")
	end
end

# ╔═╡ eb76d4e7-f3a7-425e-9c4c-300e294face3
x3

# ╔═╡ b469a288-a02b-480b-93e0-b848fdb8dc06
h3 = h(x3)

# ╔═╡ 728479a2-afc1-4dd7-9156-39bdc64da9fe
md" We are getting closer and closer to the true root at each iteration. "

# ╔═╡ a103d9a0-332d-4032-9644-de961f36dd8a
md" Let us look at a fancier version of our example using sliders... "

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
	f(x) = 0.2x^3 - 4x + 1
	
	standard_Newton(f, n, -10:0.01:10, x0, -10, 70)
end

# ╔═╡ e150a318-60ca-428d-9e9c-4c06ef37e231
md" If we have a root approximation $x_k$, we can construct a linear model of $f(x)$ using the formual for the tangent line of a differentiable function: 

$q(x) = f\left(x_{k}\right)+f^{\prime}\left(x_{k}\right)\left(x -x_{k}\right)$

Finding the root of $q(x) = 0$ is easy. We define the next approximation by the condition $q(x_{k+1}) = 0$ or stated in another way,  

$x_{k+1} = x_{k} - \frac{f(x_{k})}{f^{\prime}(x_{k})}$

Starting with an initial estimate $x_1$, this formula produces a sequence of estimates $x_2, x_3, x_4, \ldots$"

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
	
	f′(x) = ForwardDiff.derivative(f, x)   # \prime<TAB>
	
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

# ╔═╡ Cell order:
# ╟─f4226cfe-ee06-4c72-9615-fc4aedfd045c
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
# ╠═fd5842db-c6b5-44d5-b60b-eadb51913247
# ╟─3f1d4727-0f5e-4513-abb9-945eed7f9874
# ╠═d2ca4b6c-b82a-4ca7-8b7f-9274ace2c628
# ╟─7e9a56d5-b169-4cfa-afc3-276743615327
# ╠═2a1e1787-d850-4ef0-b7c7-3682d0cc8881
# ╟─62c6879b-e609-45d6-affd-c82e69356919
# ╠═91ca8fdc-d3e7-4d9a-8613-f67505805b21
# ╠═812e6f4b-c266-4e29-a5f8-465f46c1b1d2
# ╠═6a50c48b-e41c-4028-a05f-e4b001dcd786
# ╟─749d1d8a-4dd6-4509-81c0-a880d95cb8e1
# ╠═ff385d11-7916-4a68-b4d0-2ade41f252e0
# ╠═eb76d4e7-f3a7-425e-9c4c-300e294face3
# ╠═b469a288-a02b-480b-93e0-b848fdb8dc06
# ╟─728479a2-afc1-4dd7-9156-39bdc64da9fe
# ╟─a103d9a0-332d-4032-9644-de961f36dd8a
# ╟─bc9bfd83-b3b2-4f66-aaa3-0f9edf0e651a
# ╟─866cf80c-26a1-421f-98f1-702bd2de2bb4
# ╟─0af7eb48-153f-4d57-8691-82403ee3454d
# ╟─637c4017-a537-4b79-902b-a4bf508e559e
# ╠═59372abe-baf4-45f2-b7bd-339ae4dba2bb
# ╟─f560279a-1c1c-45f9-991d-ceefc575da3d
# ╟─8d02a306-7ba5-4513-a2e4-103c5637a1ac
# ╟─947f1f04-5d95-4047-8f14-1947b3178b30
# ╠═a6f17a60-2251-427b-b6b3-bacf01927ed0
# ╟─e150a318-60ca-428d-9e9c-4c06ef37e231
# ╟─3a9f8b1f-bc71-4ae1-93ba-601553d1f4bb
# ╠═5fa580cd-93fd-4a40-b816-544e7c2acf7e
# ╟─3bcbb714-2b8b-44c7-83c3-568e903a64fe
# ╠═6005e783-2f1b-42a8-9771-cd29ca264e0f
# ╠═57de80f1-b660-4589-b130-6eedd80280f7
# ╠═35bab325-8916-49b4-8f46-abd9cd77d4d9
# ╠═ceb9c187-d67b-48f2-ade0-46744eda2f4d
# ╟─ffd22348-7c3a-4832-a0e4-47bc44dd54be
# ╟─3f62345e-6359-4e2a-89f0-9680798c5f30
