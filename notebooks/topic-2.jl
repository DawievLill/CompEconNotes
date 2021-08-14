### A Pluto.jl notebook ###
# v0.15.1

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

# ╔═╡ 1c7f7f74-7c57-11eb-293a-d1be483a7ca0
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
			Pkg.PackageSpec(name="Images"), 
			Pkg.PackageSpec(name="ImageMagick"), 
			Pkg.PackageSpec(name="PlutoUI"), 
			Pkg.PackageSpec(name="HypertextLiteral"), 
			Pkg.PackageSpec(name="ForwardDiff"), 
			Pkg.PackageSpec(name="SymEngine"),
			Pkg.PackageSpec(name="Plots"), 
			Pkg.PackageSpec(name="GR"),
			Pkg.PackageSpec(name="Zygote"),
			Pkg.PackageSpec(name="FiniteDifferences"),
			])
	using Images
	using PlutoUI
	using HypertextLiteral
	using LinearAlgebra
	using ForwardDiff  # Automatic differentiation
	using SymEngine  # Symbolic differentiation
	using GR # Graphics
	using Zygote # Reverse mode automatic differentiation
	using FiniteDifferences # Finite differences
end

# ╔═╡ 7819e032-7c56-11eb-290b-23dc34edfc58
md"# Functions and Differentiation"

# ╔═╡ d88705f0-7c57-11eb-1950-bd54523e4a72
md" This session draws heavily from a course on computational thinking that is presented at MIT, which can be found [here](https://computationalthinking.mit.edu/Spring21/). Much of what we present here has been taken directly from these notes. We will start with basics on functions, in order to make sure everyone is one the same page with respect to this fundamental concept. 

Once we have explored the idea of a function representation, we move to a really cool way in which you can take derivatives. This method is called `autodiff`, which is short for automatic differentiation. With this method we can automatically compute **exact** derivatives (up to floating-point error) given only the function itself.  

This method of taking derivatives is used widely in machine learning and optimisation and has become increasingly popular over the last couple of years. An excellent reference for all things optimisation, which also uses Julia as the main code base, can be found [here](https://algorithmsbook.com/optimization/#). The main message from this session will be to invest some time into understanding the implementation of autodiff, you don't really need to fully understand how it is constructed, just how to use it."

# ╔═╡ 45aed8a2-7c59-11eb-3f69-e701041d6a30
md" ## Functions (in Julia)"

# ╔═╡ 7cfa32e4-7c58-11eb-32c0-5760739f6de4
md"""
Before we get started with automatic differentiation, let us just make sure that everyone is up to speed on the basics of functions. Remember that a function is defined as a relation which assigns to each element in the domain a **single element** in the range. A relation is a set of ordered pairs, $(x, y)$. The set of first coordinates is the domain, the set of second coordinates the range of the relation. Therefore a function is simply a mapping from values in the domain to a single value in the range. 

The first type of functions you learned about were **univariate** functions or real-valued functions functions of a single variable e.g. $f₁(x)=x^2$, $f₂(x)=\sin(x)$, $f₃(x)=x^\alpha$.

In Julia, functions can be written in short form, anonymous form, or long form.
"""

# ╔═╡ a3d72150-87e9-11eb-29de-7d2adebbf9a7
f₁(x) = x^2 # Short form

# ╔═╡ cb0ac2ce-87e9-11eb-3c26-ed50a59978be
x -> sin(x) # Anonymous form (lambda function)

# ╔═╡ ce2c915c-8b3f-11eb-0c24-5727e31119b7
md" We can give this anonymous function a name if we really wanted to." 

# ╔═╡ dabd966e-8b3f-11eb-3092-b16c6c39b80f
f₂ = x -> sin(x)

# ╔═╡ dfcad438-87e9-11eb-0fd2-218a4806a08c
# Long form

function f₃(x, α = 3)
	x^α  # Julia automatically returns the last line in the code block
end

# ╔═╡ 59934d80-9063-11eb-0956-15e084f8135a
md" Now that we have established how to think about functions, let us move on to how we can take derivatives of these functions."

# ╔═╡ dbafac10-9062-11eb-169e-cf0cbd50df46
md" ## Differentiation "

# ╔═╡ ea99e100-9062-11eb-22d6-7d295c16f30a
md" In this section we will cover different forms of differentiation. Traditionally one would do derivatives by hand and input that derivative into the computer. This is called manual differentiation. This method is slow and prone to human error. See the following example.  "

# ╔═╡ ad3c14d0-9063-11eb-153d-e3efee72a8d5
g(x) = exp(2x) - x^3

# ╔═╡ e21fa9a2-9063-11eb-389e-5f294f5762f0
g_prime(x) = 2exp(2x) - 3x^2

# ╔═╡ f618e202-9063-11eb-3f97-cb280350cbc9
md" We can evaluate the function and its derivative at certain points, but once again we want something more automatic than this. "

# ╔═╡ 2622db40-9064-11eb-156c-e9195c13dbd7
md" #### Automatic differentiation "

# ╔═╡ acf5a080-9064-11eb-0b61-55cf13f599c7
md" In this section we will quickly mention how to implement automatic differentiation, but we will not go into detail about how it actually works. The video below gives a great overview of the process of automatic differentiation. [Here](https://julia.quantecon.org/more_julia/optimization_solver_packages.html) is a more advanced set of notes on the topic of differentiable programming and automatic differentiation in Julia. "

# ╔═╡ 6a139650-9064-11eb-3f45-41506ea4e0f0
html"""

<iframe width="680" height="400" src="https://www.youtube.com/embed/wG_nF1awSSY" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

"""

# ╔═╡ 6ad1feb0-87ea-11eb-14fb-2b73c5bacf7d
md" In terms of the univariate functions we defined above it is easy to perform automatic differentiation. This method of differentiation is different from symbolic differentiaion that you will encounter in calculus or numerical differentition via first differences."

# ╔═╡ b930ded2-87ea-11eb-3f26-3d4e1597615b
md" Autodifferentiation is used in many machine learning and scientific computing applications, so it is a useful tool to know. In order to take calculate the derivative of our function $f₁(x)$ and it evaluate it at $x = 10$ we only have to do the following."

# ╔═╡ 5e60acc0-8b38-11eb-11be-3bb33c2e8c72
ForwardDiff.derivative(f₁, 10) # Derivative of x² evaluated at x = 10

# ╔═╡ cb201760-8b38-11eb-266a-0572493239ae
ForwardDiff.derivative( x -> f₃(x, 3), 10) # Derivative of x³ at x = 10

# ╔═╡ 93607132-8b38-11eb-0359-2705eb558814
md" As observed from the last equation, we could also use the anonymous function notation. In this case with the parameter set to $\alpha = 3$."

# ╔═╡ 347785f0-906f-11eb-08e9-69ba5d3e4bf3
md" One can repeat the same calculations with a reverse-mode automatic differentiation package called `Zygote.jl`. This is one of the newest and most exciting implementations of automatic differentiation in Julia and is used in the `Flux.jl` package, which is the foundational machine learning package in Julia. "

# ╔═╡ 677fe790-906e-11eb-03b4-c1545d5f7e3c
f₁'(10) # Using Zygote

# ╔═╡ fcf37cb0-906e-11eb-37f3-5bbec323e815
f₃'(10) 

# ╔═╡ c0cb6d4a-87ec-11eb-348b-e540882173e3
md" There are several packages available to perform automatic differentiation in Julia. [Here](https://juliadiff.org/) is a curated list of all actively maintained `autodiff` packages. We were using `ForwardDiff.jl` and `Zygote.jl`. You can read more about how to use `Zygote` at the [website](https://fluxml.ai/Zygote.jl/latest/). " 

# ╔═╡ 95c81368-8b3a-11eb-27bd-59e48dbc509b
md" #### Finite differences vs autodiff"

# ╔═╡ e9c48a0c-8b38-11eb-373e-3fc1ea18d52b
md" Consider the following example to see how much better the automatic differentiation representation of $f_4(x) = \sin(x)$ is compared to that of finite differences. 

Note: Remember from calculus that the derivative of $\sin(x)$ is $\cos(x)$. "

# ╔═╡ bd76ca72-8b39-11eb-0147-252776c0eddf
md" The finite difference method computes the difference between two values that differ by a finite step size. Finite differences utilises the limit definition of a derivative, namely 

$\frac{df_{4}(x)}{dx} \equiv \lim_{\epsilon \rightarrow 0}\frac{\sin(x+\epsilon) - \sin(x)}{\epsilon}$

However, we cannot use this definition directly, since we can let $t = 1/\epsilon$ and then reframe the problem such it has an infinite limit, which is shown below. 

$\frac{df_{4}(x)}{dx} \equiv \lim_{t \rightarrow \infty}\frac{\sin(x + 1/t) - \sin(x)}{1/t}$

We know from our previous session that the computer only has a finite space to store values for $t$, so to avoid overflow issues, we need to have a finite representation. A **forward difference** representation looks exactly like the formal definition but without the limit component. 

$\frac{df_{4}(x)}{dx} \approx \frac{\sin(x+\epsilon) - \sin(x)}{\epsilon}$

The slider below allows us to determine the value for magnitude of $\epsilon$. Smaller values of $\epsilon$ lead to higher precision approximation of the derivative through finite difference approximation.  

"

# ╔═╡ 75f95eb2-8b39-11eb-211f-512a656e2f36
begin
	md"""
	$(@bind e Slider(-17:-1, default=-1, show_value=true))
	"""
end

# ╔═╡ 9e029a72-8b39-11eb-0a25-6dc0aa5e1d4e
ϵ = 10.0^e

# ╔═╡ 327cf250-8b4f-11eb-16b6-e709eb78504c
md" The different approaches are finite differences, symbolic differentiation, automatic differentiation. These representations are evaluated at $x = 1$ below." 

# ╔═╡ a567131a-8b39-11eb-0769-11888a3de6b6
(sin(1+ϵ)-sin(1))/ϵ , cos(1), ForwardDiff.derivative(sin,1), sin'(1)

# ╔═╡ 4083acc0-9065-11eb-02d9-83b9411801ec
md" We see that we get different answers for the different methods. The answer also seems to depend on the value that we choose for $\epsilon$. The reason for this has to do with truncation and rounding error, topics that we mentioned in the previous session. Normally there is a tradeoff between truncation and rounding error and one would like to pick a value for $\epsilon$ that provides the lowest level of truncation and rounding error. Picking the value for $\epsilon$ can be really hard, so we can use the `FiniteDifferences.jl` package in Julia to do it for us."

# ╔═╡ 3d133d20-9070-11eb-0b1d-2533e8ab4355
forward_fdm(12, 1)(sin, 1) # 5th order forward difference

# ╔═╡ dfc3a550-9070-11eb-0ac2-21fcd1db7a0d
md" Below is a quick discussion on truncation and rounding error. This gives us an idea of why are getting our particular errors. "

# ╔═╡ 30f39ae0-9065-11eb-3acf-5766f18d038a
md" ##### Truncation error "

# ╔═╡ 9f50c068-8da8-11eb-2017-27fb9a5b87a4
md" With 'large' values of $\epsilon$ our estimate seems to be quite far off from the actual answer. This is due to truncation error. We are trying to emulate the value for $\epsilon$ going to zero, but we are choosing a nonzero value to represent $\epsilon$. Consider the following example to understand the problem.  "

# ╔═╡ 8d7d0710-9065-11eb-3b33-2be258397e78
h(x) = 5x^3

# ╔═╡ afdfaba0-9065-11eb-2dd0-d9e5380faf3c
h_prime(x) = 15x^2

# ╔═╡ c1953590-9065-11eb-2e9a-3b9cd59eb1fe
h(7)

# ╔═╡ c55aa200-9065-11eb-2fe3-dfac936d1a72
h_prime(7)

# ╔═╡ cd5f65d0-9065-11eb-391d-c3083aba15e1
h_prime_approx = (h(7 + 0.25) - h(7))/0.25

# ╔═╡ f067c042-9065-11eb-0821-719c6a92b150
TE = h_prime(7) - h_prime_approx 

# ╔═╡ 616e1e50-9067-11eb-3f55-4310dee68485
md" The question is why do we get this truncation error. Well the answer can be depicted with the following calculation. This section is optional."

# ╔═╡ 7c3c7420-9067-11eb-1aa3-236072753a44
md" We can begin with a Taylor series expansion where we have 

$f(x + \epsilon) = f(x) + \epsilon \frac{d f(x)}{dx} + \frac{\epsilon^2}{2!} \frac{d^2 f(x)}{dx^2} + \frac{\epsilon^3}{3!} \frac{d^3 f(x)}{dx^3} + \cdots$

If we use the forward difference representation to approximate our derivative then the Taylor series expansion will give the following

$\frac{f(x + \epsilon) - f(x)}{\epsilon} = \frac{d f(x)}{dx} + \left[\frac{\epsilon}{2!} \frac{d^2 f(x)}{dx^2} + \frac{\epsilon^2}{3!} \frac{d^3 f(x)}{dx^3} + \cdots \right]$

The term on the left hand side is the numerical approximation, then on the right hand side the first term is the true derivative. The component in brackets is the truncation error."

# ╔═╡ 4cc1df10-9066-11eb-0938-e56149d21210
md" ##### Rouding error "

# ╔═╡ 22188120-9065-11eb-2c2a-e9eddcfab479
md" We observe that with finite differences the answer is closer to the true value with smaller values of $\epsilon$. There is one thing that we must take into consideration here. Remember our old enemy, catastrophic cancellation!

Since we are taking differences, we can often encounter issues when the values of $\epsilon$ are really small, since this means that $\sin(x+\epsilon)$ will be close to $\sin(x)$. We need to be careful that values don't cancel out in our calculation. In fact if you run the slider back far enough you will see the cancellation occur in our example.  

Let us illustrate this point in more detail with another example. Let us calculate the derivative of $x^2$ at $x = 2$ with finite differences. "

# ╔═╡ 49de4d34-8da9-11eb-030b-59bb806f26cb
ds(ϵ, x) = ((x + ϵ)^2 - x^2)/ϵ

# ╔═╡ 78fddcce-8da9-11eb-1acd-eb95a14a2673
HTML("The deriviative with ϵ = 1e-8 is: $(ds(1e-8, 2.))")

# ╔═╡ b5691a34-8da9-11eb-09ae-af3a5f3746af
HTML("The deriviative with ϵ = 1e-12 is: $(ds(1e-12, 2.))")

# ╔═╡ b53499da-8da9-11eb-06b9-038a539e6b5a
HTML("The deriviative with ϵ = 1e-16 is: $(ds(1e-16, 2.))")

# ╔═╡ b5087744-8da9-11eb-027e-4d9008814177
HTML("The deriviative with ϵ = 1e-1 is: $(ds(1e-1, 2.))")

# ╔═╡ c932e536-8dae-11eb-1cc9-015e3ab7c8ea
((2 + ϵ)^2 - 2^2)/ϵ

# ╔═╡ b4ed0598-8da9-11eb-3dc1-d39ba25b0aef
md" So we see in this case that the derivative with $\epsilon = 10^{-16}$  actually results in subtractive cancellation and we get an answer that is completely wrong. We can solve these rounding issues by picking an optimal value for ϵ. This optimal value can be shown to be $ϵ = \max\{|x|, 1\}\sqrt{\text{eps}}$ where $\text{eps}$ refers to machine epsilon.  "

# ╔═╡ a65b9ec6-8daa-11eb-3529-1de77b221d3e
md" #### Errors and finite differences (optional) "

# ╔═╡ c3a2963a-8daa-11eb-38e5-8198e6c46ac3
md" We have to ask ourselves whether there is a way to measure the error growth rate in ϵ. In order to do this we have to perform a first order Taylor expansion of our function of interest around the value where the function is to be evaluated. In other words, Taylor expansion of $f(x)$ around $x$.

$f(x + ϵ) = f(x) + f'(x)ϵ + \mathcal{O}(ϵ^{2})$

The $\mathcal{O}(ϵ^{2})$ means that our error in the approximation grows at quadratically in ϵ, which seems strange since we took a linear approximation.
"

# ╔═╡ c3867f52-8daa-11eb-1c17-1f51a914b7a1
md" To better understand this error approximation, let us rearrange the Taylor expansion. 

$f'(x) = \frac{f(x + ϵ) - f(x)}{ϵ} + \mathcal{O}(ϵ^2)/ϵ$

This better reflects the fact that forward differences has *linearly* growing errors, since $\mathcal{O}(ϵ^2)/ϵ = \mathcal{O}(ϵ)$. This means that if we halve $ϵ$ we will also halve the error in our approximation. This ignores issues surrounding rounding though. 

"

# ╔═╡ 4ec7be62-8dad-11eb-2099-2d76bc183810
md" There are ways to improve on the accuracy of the finite difference method. One could, as an example, use central differences instead of finite differences. The central difference representation for our equation from before is:

$\frac{df_{4}(x)}{dx} \approx \frac{\sin(x+\epsilon) - \sin(x - \epsilon)}{2\epsilon}$

The error will grow much slower with this representation, as opposed to forward differences, so it seems to be more beneficial to use this. However, one needs to consider that for with central differences there are many more computations to perform, since we need to evaluate the function twice in each iteration. This means that forward differences are computationally more efficient, but less accurate. "

# ╔═╡ 3fb7e3e6-8b4d-11eb-308d-f1d31d42e184
md" #### What about symbolic differentiation? "

# ╔═╡ 1a011e5c-8b4c-11eb-1996-ed9145ec9ee7
md" An alternative to automatic differentiation and finite differences would be to use the _SymEngine.jl_ package to do symbolic differentiation. One could also use the _SymPy_ package from Python directly in Julia. For symbolic differentiation and manipulation _Mathematica_ is probably the most sophisticated, so if you have to do symbolic work, you should work in that space.  "

# ╔═╡ 57a46a0a-8b4c-11eb-1210-cf6a303f32ac
@vars s;

# ╔═╡ 5da00e52-8b4c-11eb-199c-c1b0fda27cc3
f₄ = sin(s)

# ╔═╡ 75ed19fa-8b4c-11eb-0873-4b96f4617096
diff(f₄, s)

# ╔═╡ 97c0c0a4-8b4c-11eb-218d-6fb509350f95
md" This seems to give us the answer that we were looking for. The derivative of $\sin(x)$ is $\cos(x)$."

# ╔═╡ 3e4817d0-9069-11eb-0e44-5f73ff9bfc73
md" Symbolic differentiation works quite well, but there is the problem of *expression swell*. In this case we have that derivative expressions are many times larger than the original function. You can read on the topic [here](http://www-troja.fjfi.cvut.cz/~liska/ca/node53.html). This idea is also nicely described in the video posted above. "

# ╔═╡ f603ce46-8b41-11eb-0d87-99b9c580169e
md" Next we move into the world of scalar valued multivariate functions. A scalar valued function takes one or more inputs and only retuns a single value. So this means that our function in this setting will take in multiple variables with potentially different values attached to those variables and only return a single value. A general case would be an $n$-variable scalar valued function that maps from the space $\mathbb{R}^{n}$ to $\mathbb{R}$." 

# ╔═╡ cf684c20-8b42-11eb-36e0-c318082f9f4f
md" ### Scalar valued multivariate functions "

# ╔═╡ e1cca2ee-8b42-11eb-0471-23a6523e7779
md" Let us consider an example of a scalar valued multivariate function $f_5(x): \mathbb{R}^{3} \rightarrow \mathbb{R}$,

$f_5(\textbf{x}) = 5\sin(x_1 * x_2) + 2x_2 / 4x_3$" 

# ╔═╡ 189c3176-8b44-11eb-03fe-71a833f4d5e6
md" There are multiple ways in Julia in which we can write this type of function. One can represent it as a function of many variables or a function of a vector. "

# ╔═╡ c68350a4-8b43-11eb-0113-8ff369685239
# Different ways to represent scalar valued multivariate function

begin
	f₅(x, y, z) = 5sin(x*y) + 2y/4z
	f₅(v) = 5*sin(v[1]*v[2]) + 2*v[2]/4v[3] 
end

# ╔═╡ 8711b310-8b44-11eb-3311-0b51038cab72
md" Once the code cell has been executed, you will see that a `generic function with 2 methods` has been created. In this case, depending on whether you use a vector or list with three elements, you will call a different version of $f_5$. Let us illustrate with an example. "

# ╔═╡ a3188eb2-8b44-11eb-3043-6faea3ffd81f
f₅(1,2,3) # Input is a list of three elements (x, y, z)

# ╔═╡ abde4208-8b44-11eb-1e72-29ca3336a118
f₅([1,2,3]) # Input is a vector v (with three elements)

# ╔═╡ ca5cf856-8d70-11eb-1d37-135bdd1a95c7
md" The fact that the same function name can do different things depending on the type of input is referred to as **multiple dispatch**. Most functions have this feature in Julia."

# ╔═╡ 0903d2fe-8b45-11eb-2700-292eac0f88f7
md" There is an even more efficient way of coding the portion above. Remember that programmers are lazy and don't want to copy code. You want to reuse code as much as possible. Look at the following code block and see if you can figure out how this abstraction works. I think it is quite elegant. " 

# ╔═╡ ca4cf144-8b44-11eb-11a7-9f5b1511f14f
begin
	f₆(x, y, z)  = 5sin(x*y) + 2y/4z
	f₆(v) = f₆(v[1], v[2], v[3]) # Does this part make sense to you? Think about what this is actually doing. 
end

# ╔═╡ 6685b140-8b45-11eb-08b8-6dc1fefab50b
f₆(1,2,3), f₆([1,2,3])

# ╔═╡ a83a91e8-8b4d-11eb-02f1-0506c0723d00
md" The gradient is a generalisation of the concept of derivative to multivariate functions. It provides the local slope of the function, which gives some idea of what is going to happen if take a small step in a certain direction on that function. In the single vairable setting derivatives are slopes of tangent lines, while in the multivariate case we know that the gradients point in the direction of steepest ascent of the tangent **hyperplane**. "

# ╔═╡ 7ec6d044-7c58-11eb-112d-5d8be6b4288c
md" #### Nabla! $\nabla$! "

# ╔═╡ 7ee10e48-7c58-11eb-2e3b-f5805d4af19c
md" In scientific computing and especially in machine learning one needs to take derivatives of multivariate functions in the direction of every argument. Collecting the partial derivatives in a vector provides us with the gradient. The gradient of the function $f$ evaluated at $\textbf{x}$ is normally written as $\nabla f(\textbf{x})$, where:

$\nabla f(\textbf{x}) = \left[\frac{\partial f(\textbf{x})}{\partial x_1}, \frac{\partial f(\textbf{x})}{\partial x_2}, \frac{\partial f(\textbf{x})}{\partial x_3}\right]$"

# ╔═╡ 7f14d2fa-7c58-11eb-1e88-773a06148b22
md" You can calculate this gradient explicitly in the manner specified below. This becomes a bit more difficult to do when the number of dimensions increase, but it is good to see in order to establish intuition for what we are doing. "

# ╔═╡ 7f2f1b06-7c58-11eb-038e-15bd2b4d1dbb
begin
	∂f₅∂x =  (f₅(1+ϵ, 2, 3  ) -f₅(1, 2, 3)) / ϵ
	∂f₅∂y =  (f₅(1, 2+ϵ, 3  ) -f₅(1, 2, 3)) / ϵ
	∂f₅∂z =  (f₅(1, 2,   3+ϵ) -f₅(1, 2, 3)) / ϵ
	∇f = [ ∂f₅∂x , ∂f₅∂y, ∂f₅∂z]
end

# ╔═╡ ef8bc500-8ce0-11eb-380f-618627e4d182
md" Since the method above requires many evaluations to compute, it is considered to be inefficient. A better alternative is to simply use automatic differentiation to find the gradient. In fact, most machine learning algorithms that are looking for specific loss functions are using automatic differentiation procedures."

# ╔═╡ 7efab6f4-7c58-11eb-006e-412089f351a0
ForwardDiff.gradient(f₅, [1, 2, 3])

# ╔═╡ d1606590-906e-11eb-1846-1500a526a8a8
gradient(f₅, [1, 2, 3]) # Using Zygote

# ╔═╡ f86ef080-8cdf-11eb-3540-59c2bd06a1e5
md" #### Surface plots in Julia "

# ╔═╡ d3d5e292-8cdf-11eb-1a17-8d5ceebaaf91
md" Let us draw some figures here to illustrate the function in three dimensions. We can look at the function $f_{7}(\textbf{x}) = \sin(x_1) + \cos(x_2)$. This is a relatively easy function to work with. and has a nice three dimensional representation. This part is mostly for fun, to get used to the idea of graphing in Julia. We are using the GR backend in this case, but you can use whatever plotting backend that you prefer. The plotting in Julia is quite similar to Matlab and shares some similarity with Python's `matplotlib` package. "

# ╔═╡ 7f9507f4-7c58-11eb-045d-ad36f4fb184a
x₁ = 8 .* rand(100) .- 4;

# ╔═╡ 7faea394-7c58-11eb-2529-c3881c14b364
x₂ = 8 .* rand(100) .- 4;

# ╔═╡ 50ad955c-8ce1-11eb-1e93-d7a3cc6a89fe
f₇ = sin.(x₁) + cos.(x₂); 

# ╔═╡ e1f0ef0a-8ce1-11eb-022b-adfd5162cdc2
md" Now let us plot our surface plot of the function we defined. I think it looks pretty cool! "

# ╔═╡ 88be634a-8ce1-11eb-0aa3-3f3783bc5eba
wireframe(x₁, x₂, f₇)

# ╔═╡ 6df8832c-8ce6-11eb-09bf-b15a882acc1a
md" You might have noticed in the construction of the graph that there is this `.` following the $\sin$ and $\cos$ function. You can also see this operator in previous lines of code in front of `*` and `-`. If you enact this operator you are `broadcasting` the operation across all the elements of the array."  

# ╔═╡ 4774616c-8db0-11eb-2368-cd1cc4cc884e
md" ### Broadcasting "

# ╔═╡ 5891437a-8db0-11eb-0157-8d8f9b77e1fb
md" If we take a closer look at the first line of code in the plotting section we see that we are taking the value of 8 and multiplying it with `rand(100)`. First, you need to ask yourself, what is `rand(100)` doing. Let's type it in to see. You should also consult the documentation on this command to get a better idea of what is happening when you call this function. " 

# ╔═╡ bbff790e-8db0-11eb-1c0b-4b8eff5e4b5e
rand(100)

# ╔═╡ e313133e-8db0-11eb-01eb-cd399ba999cb
typeof(rand(100))

# ╔═╡ bf8478fe-8db0-11eb-2576-45a741c4377d
md" It appears that it creates an array, in this case a **column vector**, that contains $100$ floating-point values. It chooses these values randomly from a uniform distribution in the interval $[0, 1)$."

# ╔═╡ bf9b99c6-8db0-11eb-0c14-195455bf6f0b
md" If we wanted to multiply each of the values in this vector by 8 then it seems that we could simply do scalar multiplication. In this case the scalar multiplication will work without using the `.` operator. So we can proceed with or without it. "

# ╔═╡ bfe58446-8db0-11eb-0e8d-51939d21730d
8 * rand(100)

# ╔═╡ b3f1ab8c-8db1-11eb-3cba-1792c29af5d0
8 .* rand(100)

# ╔═╡ bfff84e0-8db0-11eb-1a00-dfe087e10c2e
md" Both of these methods will work with multiplication. However, when we apply subtraction we get an error. It says that if we want element-wise subtraction we need to use broadcasting with dot syntax." 

# ╔═╡ c01b2a2e-8db0-11eb-3422-0b59816077fe
8 .* rand(100) - 4

# ╔═╡ ca24f6de-8db1-11eb-19bc-f71b8b226c39
md" This then indicates to us that the broadcasting operation allows us to subtract 4 from each of the components in the vector in an elementwise fashion. "

# ╔═╡ c03422fe-8db0-11eb-34d3-4bf88e1da37e
8 .* rand(100) .- 4

# ╔═╡ 2700f722-8db2-11eb-3993-4f7392659ef8
md" In most instances you can use the broadcast operator across all functions without worrying with it too much. Julia is quite good at figuring out when you should be broadcasting.  "

# ╔═╡ Cell order:
# ╟─1c7f7f74-7c57-11eb-293a-d1be483a7ca0
# ╟─7819e032-7c56-11eb-290b-23dc34edfc58
# ╟─d88705f0-7c57-11eb-1950-bd54523e4a72
# ╟─45aed8a2-7c59-11eb-3f69-e701041d6a30
# ╟─7cfa32e4-7c58-11eb-32c0-5760739f6de4
# ╠═a3d72150-87e9-11eb-29de-7d2adebbf9a7
# ╠═cb0ac2ce-87e9-11eb-3c26-ed50a59978be
# ╟─ce2c915c-8b3f-11eb-0c24-5727e31119b7
# ╠═dabd966e-8b3f-11eb-3092-b16c6c39b80f
# ╠═dfcad438-87e9-11eb-0fd2-218a4806a08c
# ╟─59934d80-9063-11eb-0956-15e084f8135a
# ╟─dbafac10-9062-11eb-169e-cf0cbd50df46
# ╟─ea99e100-9062-11eb-22d6-7d295c16f30a
# ╠═ad3c14d0-9063-11eb-153d-e3efee72a8d5
# ╠═e21fa9a2-9063-11eb-389e-5f294f5762f0
# ╟─f618e202-9063-11eb-3f97-cb280350cbc9
# ╟─2622db40-9064-11eb-156c-e9195c13dbd7
# ╟─acf5a080-9064-11eb-0b61-55cf13f599c7
# ╠═6a139650-9064-11eb-3f45-41506ea4e0f0
# ╟─6ad1feb0-87ea-11eb-14fb-2b73c5bacf7d
# ╟─b930ded2-87ea-11eb-3f26-3d4e1597615b
# ╠═5e60acc0-8b38-11eb-11be-3bb33c2e8c72
# ╠═cb201760-8b38-11eb-266a-0572493239ae
# ╟─93607132-8b38-11eb-0359-2705eb558814
# ╟─347785f0-906f-11eb-08e9-69ba5d3e4bf3
# ╠═677fe790-906e-11eb-03b4-c1545d5f7e3c
# ╠═fcf37cb0-906e-11eb-37f3-5bbec323e815
# ╟─c0cb6d4a-87ec-11eb-348b-e540882173e3
# ╟─95c81368-8b3a-11eb-27bd-59e48dbc509b
# ╟─e9c48a0c-8b38-11eb-373e-3fc1ea18d52b
# ╟─bd76ca72-8b39-11eb-0147-252776c0eddf
# ╟─75f95eb2-8b39-11eb-211f-512a656e2f36
# ╠═9e029a72-8b39-11eb-0a25-6dc0aa5e1d4e
# ╟─327cf250-8b4f-11eb-16b6-e709eb78504c
# ╠═a567131a-8b39-11eb-0769-11888a3de6b6
# ╟─4083acc0-9065-11eb-02d9-83b9411801ec
# ╠═3d133d20-9070-11eb-0b1d-2533e8ab4355
# ╟─dfc3a550-9070-11eb-0ac2-21fcd1db7a0d
# ╟─30f39ae0-9065-11eb-3acf-5766f18d038a
# ╟─9f50c068-8da8-11eb-2017-27fb9a5b87a4
# ╠═8d7d0710-9065-11eb-3b33-2be258397e78
# ╠═afdfaba0-9065-11eb-2dd0-d9e5380faf3c
# ╠═c1953590-9065-11eb-2e9a-3b9cd59eb1fe
# ╠═c55aa200-9065-11eb-2fe3-dfac936d1a72
# ╠═cd5f65d0-9065-11eb-391d-c3083aba15e1
# ╠═f067c042-9065-11eb-0821-719c6a92b150
# ╟─616e1e50-9067-11eb-3f55-4310dee68485
# ╟─7c3c7420-9067-11eb-1aa3-236072753a44
# ╟─4cc1df10-9066-11eb-0938-e56149d21210
# ╟─22188120-9065-11eb-2c2a-e9eddcfab479
# ╠═49de4d34-8da9-11eb-030b-59bb806f26cb
# ╟─78fddcce-8da9-11eb-1acd-eb95a14a2673
# ╟─b5691a34-8da9-11eb-09ae-af3a5f3746af
# ╟─b53499da-8da9-11eb-06b9-038a539e6b5a
# ╟─b5087744-8da9-11eb-027e-4d9008814177
# ╠═c932e536-8dae-11eb-1cc9-015e3ab7c8ea
# ╟─b4ed0598-8da9-11eb-3dc1-d39ba25b0aef
# ╟─a65b9ec6-8daa-11eb-3529-1de77b221d3e
# ╟─c3a2963a-8daa-11eb-38e5-8198e6c46ac3
# ╟─c3867f52-8daa-11eb-1c17-1f51a914b7a1
# ╟─4ec7be62-8dad-11eb-2099-2d76bc183810
# ╟─3fb7e3e6-8b4d-11eb-308d-f1d31d42e184
# ╟─1a011e5c-8b4c-11eb-1996-ed9145ec9ee7
# ╠═57a46a0a-8b4c-11eb-1210-cf6a303f32ac
# ╠═5da00e52-8b4c-11eb-199c-c1b0fda27cc3
# ╠═75ed19fa-8b4c-11eb-0873-4b96f4617096
# ╟─97c0c0a4-8b4c-11eb-218d-6fb509350f95
# ╟─3e4817d0-9069-11eb-0e44-5f73ff9bfc73
# ╟─f603ce46-8b41-11eb-0d87-99b9c580169e
# ╟─cf684c20-8b42-11eb-36e0-c318082f9f4f
# ╟─e1cca2ee-8b42-11eb-0471-23a6523e7779
# ╟─189c3176-8b44-11eb-03fe-71a833f4d5e6
# ╠═c68350a4-8b43-11eb-0113-8ff369685239
# ╟─8711b310-8b44-11eb-3311-0b51038cab72
# ╠═a3188eb2-8b44-11eb-3043-6faea3ffd81f
# ╠═abde4208-8b44-11eb-1e72-29ca3336a118
# ╟─ca5cf856-8d70-11eb-1d37-135bdd1a95c7
# ╟─0903d2fe-8b45-11eb-2700-292eac0f88f7
# ╠═ca4cf144-8b44-11eb-11a7-9f5b1511f14f
# ╠═6685b140-8b45-11eb-08b8-6dc1fefab50b
# ╟─a83a91e8-8b4d-11eb-02f1-0506c0723d00
# ╟─7ec6d044-7c58-11eb-112d-5d8be6b4288c
# ╟─7ee10e48-7c58-11eb-2e3b-f5805d4af19c
# ╟─7f14d2fa-7c58-11eb-1e88-773a06148b22
# ╠═7f2f1b06-7c58-11eb-038e-15bd2b4d1dbb
# ╟─ef8bc500-8ce0-11eb-380f-618627e4d182
# ╠═7efab6f4-7c58-11eb-006e-412089f351a0
# ╠═d1606590-906e-11eb-1846-1500a526a8a8
# ╟─f86ef080-8cdf-11eb-3540-59c2bd06a1e5
# ╟─d3d5e292-8cdf-11eb-1a17-8d5ceebaaf91
# ╠═7f9507f4-7c58-11eb-045d-ad36f4fb184a
# ╠═7faea394-7c58-11eb-2529-c3881c14b364
# ╠═50ad955c-8ce1-11eb-1e93-d7a3cc6a89fe
# ╟─e1f0ef0a-8ce1-11eb-022b-adfd5162cdc2
# ╠═88be634a-8ce1-11eb-0aa3-3f3783bc5eba
# ╟─6df8832c-8ce6-11eb-09bf-b15a882acc1a
# ╟─4774616c-8db0-11eb-2368-cd1cc4cc884e
# ╟─5891437a-8db0-11eb-0157-8d8f9b77e1fb
# ╠═bbff790e-8db0-11eb-1c0b-4b8eff5e4b5e
# ╠═e313133e-8db0-11eb-01eb-cd399ba999cb
# ╟─bf8478fe-8db0-11eb-2576-45a741c4377d
# ╟─bf9b99c6-8db0-11eb-0c14-195455bf6f0b
# ╠═bfe58446-8db0-11eb-0e8d-51939d21730d
# ╠═b3f1ab8c-8db1-11eb-3cba-1792c29af5d0
# ╟─bfff84e0-8db0-11eb-1a00-dfe087e10c2e
# ╠═c01b2a2e-8db0-11eb-3422-0b59816077fe
# ╟─ca24f6de-8db1-11eb-19bc-f71b8b226c39
# ╠═c03422fe-8db0-11eb-34d3-4bf88e1da37e
# ╟─2700f722-8db2-11eb-3993-4f7392659ef8
