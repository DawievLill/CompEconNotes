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

# ╔═╡ 7fffd413-6ebe-49d2-b426-1428aee5ae26
begin
	using ApproxFun
	using Colors
	using Images
	using Interpolations
	using LaTeXStrings
	using NLsolve
	using Optim
	using Parameters
	using Plots
	using PlutoUI
	using QuantEcon
	using Random
end

# ╔═╡ c428b2b9-a4f4-4ddc-be20-d25cb14c9cf7
html"""
<style>
  main {
    max-width: 900px;
  }
</style>
"""

# ╔═╡ 56121f90-e33f-11eb-1211-8578ae4eb05d
md" # Dynamic programming I"

# ╔═╡ 058c997b-390b-4b8e-b203-20879e0673a4
md" ## Still quite incomplete and incoherent. Needs work! "

# ╔═╡ f6852b65-d493-4e2c-a8b6-517993c93ac8
md" The primary sources for this session is the numerical methods course by Florian Oswald and the set of notes at QuantEcon. Dynamic programming is a method to solve dynamic problems in economics, and other disciplines. It is quite sseful in several areas of economics, but I have mostly used it in macroeconomic applications. It is therefore natural to me that most of the examples illustrated here will be from a macroeconomic context, but I am open to suggestions for problems in labour economics, microeconomics, etc. We will start with the simplest dynamic programming problem, referred to as the shortest path problem and then discuss the cake eating problem and optimal growth problem.

In this session we cover some basic discrete deterministic dynamic programming problems. We will first look at models with a finite time horizon and then move on to models with an inifinite time horizon. The techniques we will utilise include backward induction and two versions of value function iteration. "

# ╔═╡ 3e79f184-a3e8-44a0-9609-19827522d86a
md" ## Shortest path problem "

# ╔═╡ 340698bc-e27a-442b-a6fe-4b7c7ed3bd16
md"""
Let's start by looking at the following problem. Create a random matrix and follow paths on it. The paths start at one of the square on the top, and can only go downwards, either South-East, South, or South-West.

We will add up the numbers visited along each path. Our goal is to find the path that has the smallest sum. So this is indeed an optimization problem: we want to **minimize** the sum along these particular paths.
"""

# ╔═╡ 08bf10d7-0221-40e8-99ac-4087ad459add
md"""
n = $(@bind n Slider(2:12, show_value = true, default=8))
"""

# ╔═╡ 6071a006-d300-487b-835e-bc4b8487d458
M = rand( 0:9, n, n);

# ╔═╡ 621779a2-fcbd-4e55-8b36-447733b92ea5
begin
	struct Paths
	    m::Int
	    n::Int
	end
	
	Base.iterate(p::Paths) = fill(1,p.m), fill(1,p.m) #start the iteration with 1's
	
	Base.IteratorSize(::Type{Paths}) = SizeUnknown()
	
	function Base.iterate(p::Paths, state)
		if state ≠ fill(p.n,p.m) # end when each row has an n
	      newstate = next(state,p.n)
	      return newstate, newstate
	    end
	end
	
	
	function next(path,n)
	    k = length(path)
		# start from the end and find the first element that can be updated by adding 1
	    while  k≥2 && ( path[k]==n || path[k]+1 > path[k-1]+1 )
	        k -= 1
	    end   
	    path[k] +=1 #add the one then reset the following elements
	    for j = k+1 : length(path)
	        path[j] = max(path[j-1]-1,1)
	    end
	    return(path)
	end
	

	
	function allpaths(m,n)
     v=Vector{Int}[]
	 paths = Paths(m,n)
     for p ∈ paths
        push!(v,copy(p))
    end
    v
	end
end

# ╔═╡ eea67d06-ac30-4e49-b307-b6b4691f1d76
begin
	paths = allpaths(n,n)
	numpaths = length(paths)
	md"There are $numpaths paths to check for this example"
end

# ╔═╡ f6c5bc76-b0d8-46e0-9994-b6163eee9704
md"""
Path $( @bind whichpath Slider(1:numpaths, show_value=true) )
"""

# ╔═╡ fd77c87a-ac60-4507-a251-19cb4b82088a
begin
	winnernum = argmin([sum( M[i,p[i]] for i=1:n) for p∈paths]);
	winner = paths[winnernum];
	winnertotal = sum( M[i,winner[i]] for i=1:n);
	md" The winner total is $winnertotal and the winner number is $winnernum"
end

# ╔═╡ 0c7a94ba-dd82-4015-9c39-adaec2e424f5
let
	
	path = paths[whichpath]
	values = [ M[i,path[i]] for i=1:n]
	nv = length(values)
	thetitle = join([" $(values[i]) +" for i=1:nv-1 ]) * " $(values[end]) = $(sum(values))";
	
	
	rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
	plot()
	for i=1:n, j=1:n
	   plot!(rectangle(1,1,i,j), opacity=.2, color=[:red,:white][1+rem(i+j,2) ])
	   
	end
	for i=1:n, j=1:n
	  annotate!((j+.5),n+2-(i+.5), M[i,j])
	end
	
	# The winning path 
		for i = 1:n-1
		plot!([ winner[i+1]+.5, winner[i]+.5  ],[n-i+.5, n-i+1.5], color=RGB(1,.6,.6),  linewidth=4)
	end
	
	
	for i = 1:n-1
		plot!([ path[i+1]+.5, path[i]+.5  ],[n-i+.5, n-i+1.5], color=:black,  linewidth=4)
	end
	
	plot!(xlabel="winner total = $winnertotal", xguidefontcolor=RGB(1,.5,.5))

	
	for i=1:n,j=1:n
		plot!(rectangle(.4,.4,i+.3,j+.3), opacity=1, color=RGB(0,1,0), linewidth=0,fillcolor=[RGBA(1,.85,.85,.2),:white][1+rem(i+j,2)])
	end
	plot!(title=thetitle)
	plot!(legend=false, aspectratio=1, xlims=(1,n+1), ylims=(1,n+1), axis=nothing)
end

# ╔═╡ d2de514f-ae9d-4d89-ba5e-e5a2d4c55081
md" This example is useful to get an idea of what dynamic programming is about. Now let us move to more complex problems. We will also need to talk about some of the theory surrounding dynamic programming. We will not delve too deeply into the theory, but rather provide sources that give a good overview of theorethical constructs. Dynamic programming can become quite technical and it is isnt witin the scope of these sessions to provide a rigorous training in these methods. "

# ╔═╡ 1fa1bdc2-a076-42be-9540-1d20ebcffeb5
md" ## Discrete dynamic programming "

# ╔═╡ e5e5bf22-a1d3-4bc3-8968-aa1d0545c2c8
md" This class of problems is more closely related to economics. It is a consumption-savings problem and is usefull both in micro and macroeconomics. The notes presented here are almost an exact copy of the ones used by **Florian Oswald**. "

# ╔═╡ b844d7eb-b152-4c37-9922-a0223f3b5589
md" ### Finite time cake eating problem [🍰 ≅ 😁]"

# ╔═╡ ba5f4d25-2ed0-4d5b-8f84-25f10dc0ad82
md" We want to maximise total utility given some constraint. Our first intuition is to solve this as a constrained optimisation problem. This is entirely plausible. However, dynamic programming offers some distinct advantages, especially when it comes to numerical methods. "

# ╔═╡ e5ef57d1-e00b-4005-913a-7d6579c9c59c
md" In this cake eating (🍰) problem you are given a certain budget of resources $K$ to spend (consume) over $t$ dates. In this example $K$ refers to the size of the cake. You take an action $c_t \in \mathbb{N}$ at a specific date $t$. The action refers to how much to spend (or cake to eat) in the specified period. We let $U_t(c_t)$ be the reward / utility that results from the action that you take in period $t$.

In this example the action is the consumption of cake, so this can be viewed as a consumption-savings problem. However, to keep things general we will refer to action instead of consumption. Now let us frame the problem more formally. "

# ╔═╡ 398a316b-fc53-444a-aea7-169231c95f4b
md" #### Budgeting problem "

# ╔═╡ b0ca2d8d-dfa9-4fed-b34d-90167e950153
md" The goal is to maximise total utility subject to a resource constraint.

$$\max _{c_t} \sum_{t \in \mathcal{T}} \beta^{t} U_{t}\left(c_{t}\right) \quad \text{s.t.} \quad \sum_{t \in \mathcal{T}} c_{t} = K$$ 

For this example we will not allow borrowing, so $c_t \geq 0$. For this problem our control variable is the action $c_t$. We want to choose the best action in each period. To choose the best action available in each period we want a list of all the values for $c_t$ in each period given the resources available in that period. We seek a best policy that prescribes the action that should be taken at each state at each point in time to maximise rewards. We can represent the values with a function, called the **value function**.

This value function is a function of the resources available in that period. In more mathematical terms, 

$V_{t}(K_{t}) = \max _{c_{1}, c_{2}, \ldots, c_{T}}\left\{U_{t}(c_{1})+\beta U_{t}(c_{2})+\beta^{2} U_{t}(c_{3})+\cdots+\beta^{T} U_{t}(c_{T})\right\}$

This value function specifies the maximium attainable sum of current and furure rewards given the state and time period. The value of resources evolve according to the following equation

$K_{t+1} = K_t - c_t = g(K_t, c_t)$

which is known as the **transition function**. While $c_t$ is referred to as the control variable, $K_t$ represents the state variable. The state and action spaces are considered finite in this example. Since our problem is deterministic in nature, the next period's state is known with certainty once the current period's state and action are known (which is reflected in the transition function).  "



# ╔═╡ db4c8010-a33f-4e78-bea4-b3168c5ce72b
md" #### Bellman equation "

# ╔═╡ 246e0451-111d-423c-a76f-eb4339b2d285
md" Dynamic programming is an analytical approach in which a multiperiod model is decomposed into a sequence of two period models. The idea of dynamic programming is based on the Principle of Optimality,

> An optimal policy has the property that, whatever the initial state and decision are, the remaining decisions must constitute an optimal policy with regard to the state resulting from the first decision.

Bellman's Principle implies that the value function must satisfy Bellman's recursion equation. This relationship between $V_{t}(K_{t})$ and $V_{t+1}(K_{t+1})$ is given by 

$$V_{t}\left(K_{t}\right)=\max _{0 \leq c_{t} \leq K_{t}}\left[U_{t}\left(c_{t}\right)+\beta V_{t+1}\left(g\left(K_{t}, c_{t}\right)\right)\right]$$ 

The value of having $K_t$ resources left in period $t$ is the value of optimizing current spending (the $\max U_t(c_t)$ part) plus the value of then having $g(K_t, c_t)$ units left going forward (i.e. the next state).

Practically, given $K_t$ we could just try out different $c_t$ and see which gives the highest value. How do we determine $V_{t+1}(\cdot)$? 🤔"

# ╔═╡ 728e0ddc-f482-4d9a-b6f2-41f7cf6b8619
md" #### Backward induction "

# ╔═╡ 3f18f963-ad85-4c75-ab8c-64ac3cbf312c
md" One way in which we could do this in our **finite time setup** is to specify that time stops at $t = T$ and then use a technique called backward induction. This method is available in all finite time horizon problems that allow dynamic programming structure. The value of leaving resources would then be zero, i.e. $V_{T+1}(K) = 0$. It is not optimal to leave cake on the table! Therefore the last period is going to be

$$V_{T}\left(K_{T}\right)=\max _{0 \leq c_{T} \leq K_{T}} U_{T}\left(c_{T}\right)$$

To maximise this we simply choose the value for $c_T$ that provides the highest payoff. This period's maximisation problem is then quite easy to solve. Now consider the problem at $T-1$: 

$V_{T-1}\left(K_{T-1}\right)=\max _{0 \leq c_{T-1} \leq K_{T-1}}\left[U_{T-1}\left(c_{T-1}\right)+\beta V_{T}\left(g\left(K_{T-1}, c_{T-1}\right)\right)\right]$

This is also easy to solve. We know the answer to $V_{T}\left(K_{T}\right)$ from the previous iteration. This means that we know the future. For period $T-1$ we are now left with solving the final piece of the puzzle, namely 

$\max _{0 \leq c_{T-1} \leq K_{T-1}} U_{T-1}\left(c_{T-1}\right)$

This iterative procedure will work all the way till we reach the first period. If the control variable is discrete we do not need assumptions of the payoff function (other than it is finite).

"

# ╔═╡ fc6640d2-f2de-4d09-8787-39954cc835b4
md" ### Backward induction example ($T = 2$)"

# ╔═╡ a897d516-60e9-4b96-9f1f-f7dea164f392
md" Below we provide an example to illustrate the idea of bakward induction. Imagine a case with two periods where we have that the per period utility function is given by $U_t(c_t) = \sqrt{c_t}$. 

Optimal choice for the different periods are as follows, 

Period $2$: $V_2(K_2) = \sqrt{K_2}$ because $c_2(K_2)=K_2$ where $K_2$ are the available resources in period $2$.

Period $1$: $V_1(K_1) = \max_{c_1} \sqrt{K_1} + \beta V_2(K_1 - c_1)$"

# ╔═╡ f17bd25f-a459-44b6-862e-5d8d7b5062c9
md" Let us look at things from a numerical perspective. First, we construct a grid over period $2$ resources. "

# ╔═╡ be62b6ff-a947-47ab-8905-6ac631abf8e6
@bind max_K Slider(0:100, show_value = true, default = 5)

# ╔═╡ cd9f3fd5-307a-4e1d-b92f-1bf7d441d483
md"  In this case, we have a maximum of $max_K resources in period 2. "

# ╔═╡ d9a659a2-d0dd-45bb-bb5b-fb85f0208819
grid_K2 = range(0, stop = max_K, step = 1) |> collect # Using a pipe operator (similar to one from R)

# ╔═╡ a25b4e49-46c0-4c7d-b80e-e108614ec288
C2 = grid_K2 # Consume all available resources at each grid point

# ╔═╡ 11dc1b7f-866e-42e6-b4ef-0b9d3b34366f
V2 = sqrt.(C2) # Apply the explicit function form for utility function

# ╔═╡ f2075369-02b7-4262-a543-91da85169e65
plot(grid_K2, V2, xlab = "Resources (Cake size)", ylab = "Value",label = L"V_2", m = (:circle), leg = :bottomright, line = 1.5)

# ╔═╡ 8b8c03cb-b8c5-45a8-9d45-a6730f0cf420
md" Say that we want to find the value of $K_2 = 2$. We have to reference the relevant index point for 2. The index in $K_2$ is 3 since we have 1-based indexing in _Julia_. Next we find the optimal consumption in period $1$ given a level of $K_1$."

# ╔═╡ ed53fb3c-9fe6-47d5-b465-3af3369f9177
@bind K1 Slider(1:max_K, show_value = true, default = 3) # Define different values for the level of cake in period 1.  

# ╔═╡ d609467a-6898-430f-9fc7-70cecf9cb3d0
C1 = 0 # Is this value optimal? Should not be, in fact it would be our least likely consumption value. 

# ╔═╡ be24ede3-87c7-402c-adbc-2f003e5099a0
β = 0.99 # Discount parameter

# ╔═╡ 8e14ac14-af81-4bd6-bdac-f850bb9020e9
V1_guess = fill(NaN, K1 + 1) # Initialise with NaN. Could also do this with zeros -- zeros(K1 + 1)

# ╔═╡ 9300d7d9-9563-4f85-b835-54a8db884317
K2 = K1 - C1 # How much cake is left in the second period after we made our decision in the first period. 

# ╔═╡ 53bdea85-8d59-4669-85ca-04119131064a
begin
	V1_guess[1] = sqrt.(C1) .+ β .* V2[K2 + 1] # Fill the first value of our value function with the Bellman equation
	V1_guess
end

# ╔═╡ eff5cdcf-3a07-4d08-9152-39c0b71eac3b
md" The function below encapsulates the logic that we have explored above in a neat container. _Julia_ is mostly meant for functional programming, so try to write with functions in mind. "

# ╔═╡ b6cbcf47-669d-4ae3-aa79-b83d59e04d9d
function backward_induction1(K1, max_K, grid_K2 = range(0, max_K, step=1) |> collect,
		V2 = sqrt.(grid_K2); β = 0.99)
	
	# Function not written to be elegant, mostly for exposition and to be self-contained. 
	
	V1_guess = zeros(K1 + 1) # Fill with (K1 + 1) zeros
	
	for i_guess in 0:K1 # Running this loop over the state space
		C1 = i_guess # Start with a guess of zero consumption and then moving on to higher levels.
		
		# NB In this routine the value for K2 is going to be the same as its position in the index
		K2 = K1 - i_guess # Cake left after consumption in the first period gets consumed in second period. 
		V1_guess[i_guess + 1] = sqrt.(C1) .+ β .* V2[K2 + 1]  # Value associated with cake consumed across both periods (taking control variables C1 and K2 into account)
	end
	
	idx_max = argmax(V1_guess) # Pick out the position of the argument that maximises the value function
	# Note, max function provides the actual value, while argmax gives index value. 
	V1, C1 = round(V1_guess[idx_max], digits = 3), (idx_max - 1)
	#md" The optimal consumption is $C_1 \rightarrow$ $C1 with associated value of $V1"
	#plot(V1_guess, xlab = "Consumption (Optimal at $C1)", ylab = "Value (Optimal at $V1)", label = L"V_1", line = 2)
end

# ╔═╡ 437462d4-2dd1-4ef9-aa7d-b0cb64bef455
md" Below we provide a graph of the different values associated with consumption for different levels of cake in period 1. "

# ╔═╡ 18b69513-700f-4300-b578-4742112a23ca
@bind K1_new Slider(0:max_K, show_value = true, default = 3)

# ╔═╡ 53c21b1b-088f-4604-9246-fc6796c3e256
begin
	v_collect = [backward_induction1(i, max_K)[1] for i in 0:K1_new] 
	plot(v_collect, xlab = L"K_1", ylab = "Value",label = L"V_1", m = (:circle), leg = :bottomright, line = 2)
	
end

# ╔═╡ fd4d1601-5a71-4ce1-a3bd-5d34776502f9
# @bind highR Slider(2:200,show_value = true, default = 20)

# ╔═╡ c8533761-ae1d-47a2-9860-ad88ecca2d34
md" ### Backward induction example ($T > 2$)"

# ╔═╡ 547def39-8a65-48bc-bb30-d28d016cd6fa
md" Let us consider the same example but with an arbitrary number of periods (even approaching infinity if we wish). We will start with 500 time periods. The solution remains largely the same."

# ╔═╡ 7b63094f-ddde-4dff-911a-58d8fc7226a2
begin
	# final period T
	points = 500; # The number of points on the grid (let us keep this fixed for now)
	
	# Lower and upper limits of the state space. In this case this represents the size of the cake. 	
	lowK = 1e-6; # Lowest possible of value of K (don't select 0, their are numerical consequences)
	highK = 30.0 # Highest value that K can take, given that consumption is zero. 
	
	# Log and then exponentiate for more points towards zero, which makes a nicer plot
	Kspace = exp.(range(log(lowK), stop = log(highK), length = points)); # The state space
	CT = Kspace; # Consume whatever is left -- Consume all resources in this case, since it is all left.
	VT = sqrt.(CT);  # Utility of that consumption -- This is known for VT, since there is no VT + 1
end

# ╔═╡ c5ea0324-8e8c-49d1-88f6-3abc4278aefa
function plotVT()
	plot(Kspace, VT, xlab = "K", ylab = "Value",label = L"V_T", m = (:circle), leg = :bottomright)
end

# ╔═╡ 4716e611-3342-4c3e-afc2-5afd8de2fc15
plotVT()

# ╔═╡ aa534428-e995-4885-bfeb-135ad92129a4
md" As we stated in our discussion, now look to an answer for $V_{T-1}$, given that we have the function $V_{T}$."

# ╔═╡ b40713e0-54db-4b49-8d35-a5b68f71ac89
begin
	# solving problem for period T-1 
	# now for each value in the state space [Kspace], we need to decide how much to consume
	
	# Initialise some vectors 
	w = zeros(points) # temporary vector for each choice of K -- could also use NaN as before. 
	VT_1 = zeros(points) # optimal value in T-1 at each state of K 
	ix = 0 # optimal index of action in T-1 at each state of K (important to include!)
	CT_1 = zeros(points) # optimal action in T-1 at each state of K

	for (ik, k) in enumerate(Kspace) # for all possible K-values (enumerate loops over entire state space and also creates an index ik) -- index component is NB!
        
		# in our inner loop we want to loop over all possible action choices
        for (ic, c_choice) in enumerate(Kspace)
            if k <= c_choice   # check whether that choice is feasible 
                w[ik] = -Inf   # If the choice is feasible, then assign this value to first element
            else
				k1 = k - c_choice  # tomorrow's K 
				jk = argmin(abs.(Kspace .- k1))  # index of that value in Kspace
                w[ic] = sqrt(c_choice) + VT[jk]   # value of that c_choice
            end
        end
        # find best action
        VT_1[ik], ix = findmax(w) # stores value and policy (index of optimal choice)
		CT_1[ik] = Kspace[ix]  # record optimal action level
		
    end
end

# ╔═╡ 7854bba6-88c3-4f19-b7b1-1ae8487ac9bc
let
	pv = plotVT()
	plot!(pv, Kspace, VT_1, label = L"V_{T-1}", m = (:star))
end

# ╔═╡ b71abd7b-191f-4044-91a9-d0f4e4883f88
plotaT() = plot(Kspace, CT ,label = L"a_T",leg = :topleft,ylab = "action",xlab = "K", line = 2)

# ╔═╡ 6f2eac77-4ad3-4da7-95a7-c958ccf270cd
let
	pa = plotaT()
	plot!(pa, Kspace, CT_1, label = L"a_{T-1}", line = 2)
end

# ╔═╡ d81d8f8a-ac8a-4d01-89c8-bb62960a69bc
# period T-1 till 1
# each period is the same now!
# that calls for a function!
"""
	Vperiod(grid::Vector,vplus::Vector)

Given a grid and a next period value function `vplus`,
calculate current period optimal value and actions.
"""
function Vperiod(grid::Vector,vplus::Vector)
	points = length(grid)
	w = zeros(points) # temporary vector for each choice or R'
	Vt = zeros(points) # optimal value in T-1 at each state of R
	ix = 0 # optimal action index in T-1 at each state of R
	at = zeros(points) # optimal action in T-1 at each state of R
	
	for (ir,r) in enumerate(grid) # for all possible R-values
		# loop over all possible action choices
		for (ia,achoice) in enumerate(grid)
			if r <= achoice   # check whether that choice is feasible
				w[ia] = -Inf
			else
				r1 = r - achoice  # tomorrow's R
				jr = argmin(abs.(grid .- r1))  # index of that value in Rspace
				w[ia] = sqrt(achoice) + vplus[jr]   # value of that achoice
			end
		end
		# find best action
		Vt[ir], ix = findmax(w) # stores Value und policy (index of optimal choice)
		at[ir] = grid[ix]  # record optimal action level
	end
	return (Vt, at)
end

# ╔═╡ 7529e1a7-c41c-49ec-8885-6c7c14093700
function backwards(grid, nperiods)
	points = length(grid)
	V = zeros(nperiods,points)
	a = zeros(nperiods,points)
	V[end,:] = sqrt.(grid)  # from before: final period
	a[end,:] = collect(grid)

	for it in (nperiods-1):-1:1
		x = Vperiod(grid, V[it+1,:])	
		V[it,:] = x[1]
		a[it,:] = x[2]
	end
	return (V,a)
end

# ╔═╡ b22748fb-9acb-4105-9fbe-654daf34dbf4
@bind nperiods Slider(2:20,show_value = true, default = 5)

# ╔═╡ 86f2d8ca-ea32-4558-8133-bce4789ad105
V,C = backwards(Kspace, nperiods);

# ╔═╡ b1e0e8e1-858f-4f12-a272-7197a117af25
let
	cg = cgrad(:viridis)
    cols = cg[range(0.0,stop=1.0,length = nperiods)]
	pv = plot(Kspace, V[1,:], xlab = "K", ylab = "Value",label = L"V_1",leg = :bottomright, color = cols[1])
	for it in 2:nperiods
		plot!(pv, Kspace, V[it,:], label = L"V_{%$(it)}", color = cols[it])
	end
	pv
end

# ╔═╡ 23bdbe50-8544-4bb9-9e91-ba7397ca4db2
let
	cg = cgrad(:viridis)
    cols = cg[range(0.0,stop=1.0,length = nperiods)]
	pa = plot(Kspace, C[1,:], xlab = "K", ylab = "Action",label = L"c_1",leg = :topleft, color = cols[1])
	for it in 2:nperiods
		plot!(pa, Kspace, C[it,:], label = L"c_{%$(it)}", color = cols[it])
	end
	pv = plot(Kspace, V[1,:], xlab = "K", ylab = "Value",label = L"V_1",leg = :bottomright, color = cols[1])
	for it in 2:nperiods
		plot!(pv, Kspace, V[it,:], label = L"V_{%$(it)}", color = cols[it])
	end
	plot(pv,pa, layout = (1,2))
end

# ╔═╡ 8ec82051-79b0-4cd7-8585-163ffde2b290
bar(1:nperiods,C[:, end], leg = false, title = "Given K_t = $(Kspace[end]), take action...",xlab = "period")

# ╔═╡ 616d4193-6ba7-450e-ad3e-d6b07cb90d2e
md" ## Optimal growth problem" 

# ╔═╡ bfeff134-00d8-4ca5-9902-6b8dde5facd3
md" Another example where dynamic programming is often used is the optimal growth problem. This is a simple extension of the cake eating problem that we encountered in the previous section. I will not go into detail on the dynamic programming theory, we will simply take a look at the model setup and how to solve the given problem. " 

# ╔═╡ f4654f14-6bef-452f-a9af-ec62613ad20b
md" Payoffs over time for the optimal growth model are

$$\sum_{t=1}^{\infty}\beta^{t}U\left(K_{t},c_{t}\right)$$

where $\beta<1$ is a discount factor, $K_{t}$ is the state, $c_{t}$ is the control. As we stated in the previous section, the state evolves according to $K_{t+1} = g(K_t, c_t)$. History of decisions is contained in the state variable. For this example we assume that both the payoff $u$ and transition function $g$ are stationary (do not depend on time). Problem can then be written as 

$$V(K)=\max_{K'\in\Gamma(K)}U(K,K')+\beta v(K')$$

In this case $\Gamma(K)$ is the constraint set (or feasible set) for $K'$ when the current state is $K$. The setup for the deterministic optimal growth model is then as follows, 

$$\begin{aligned}
   V(K) &= \max_{0<K'<f(K)} U(f(K) - K') + \beta V(K')\\
  f(K)  & = K^\alpha
\end{aligned}$$

where $K_0$ is given. The numerical approximation that we are going to implement is the following:

$$V(K_i) = \max_{i'=1,2,\dots,n} U(f(K_i) - K_{i'}) + \beta V(K_{i'})$$"

# ╔═╡ 5a3df38d-2973-421e-8da0-a190653af20f
md" ### Value funtion iteration "

# ╔═╡ e8e3eb80-0c0f-4a08-b963-eb554cac8f0b
md" The first method that we will use to solve this problem is value function iteration. We have already used something similar to this method for the cake eating problem, but let illustrate it for the optimal growth problem as well. For value function iteration we want to find the fixed point of the functional equation by iteration. 

We iterate untill the distance (measured in whatever way we specify) between consecutive iterations becomes small. The notion of smallness is also arbitrarily defined. The reason that we can implement this method is as a result of the contraction mapping theorem and its application to the Bellman operator. If you are interested in the theoretical motivation, I would advise you start by reading the relevant material on the QuantEcon website and also look at textbook treatments such as Lucas and Stokey. "

# ╔═╡ 8f7abbcf-61c5-4c4f-92f4-41c02a6aaa41
md" We start with discrete dynamic programming, which means that we will look at a discrete state and action space. We solve the functional problem on the real line (which has continuous support) on a finite set of grid points only. This method is quite robust and easy to understand, but is known for being slow. Eventually, given an increased number of grid points to evaluate across and enough iterations one should approximate the true solution, but this can take a lot of time. This is especially true if the state space is mutlidimensional, which brings in the curse of dimensionality. " 

# ╔═╡ 51c0c398-81d0-4e93-9907-34f3f39be2e4
md" While the general steps provided above explain the general idea quite well, it is not practical as an algorithm to solve our problem. Let us expand more on these steps to provide us with an algorithm that we can take to the computer. The basic outline of the VFI algorithm that we need to follow is the following. 

1. Set the appropriate parameter values for the model
2. Define a grid of admissable values for the state variable $K$
3. Provide initial guess for the value funtion $V$ and choose stopping criterion $\varepsilon > 0$
4. Start the iteration, for each $K_i$, $i = 1, \ldots, N$ compute

$$W(K_i) = \max_{0 \leq K^{\prime} \leq f(K)} \{U(f(K_{i}) - K^{\prime}) + \beta V(K^{\prime})\}$$

5. Stop if $d(W, V) < \varepsilon$, otherwise set $V = W$ and go back to step $4$
6. Plot value and policy functions
7. Compare with analytical solutions (provided above) "

# ╔═╡ b2f96156-bd23-4767-8ec8-38215a8f4b63
md" #### VFI implementation "

# ╔═╡ 60912ad4-c766-4b70-b843-baed710e29ce
begin
	# Step 1: Set up the appropriate parameter values for the model
	α               = 0.65
	#β 				= 0.96 		# Discount factor (already defined in global scope of this notebook)
	γ 				= 1.5 		# Degree of relative risk aversion
	maxit           = 300 		# Maximum number of iterations
end

# ╔═╡ adc503f5-6f65-4f6a-ba10-965f78f52e9e
begin
	# Step 2: Define a grid of values for the state variable 
	K_grid_min 		= 1e-3 		# Don't include zero here. Numerically more stable 
	K_grid_max 		= 2.5 		# Size of the cake
	K_grid_size 	= 120 		# Size of the grid
	
	# Set up the grid
	K_grid          = range(K_grid_min, K_grid_max, length = K_grid_size)
end

# ╔═╡ cd71f975-5058-46a4-80d8-0d410f19cbf1
begin 
	# Step 3: Initial guess for the value function and stopping criterion
	n₂ = K_grid_size  			# Number of grid points
	v₀ = zeros(n₂) 				# Initial guess is an array of zeros that have the same dimension as number of grid points.  
	
	# Tolerance level for stopping criterion
	ε = 1e-9
end

# ╔═╡ a8cbe356-218d-47f1-b62a-8fea4f6c1835
md" Quick word about the Bellman operator and contraction mapping theorem here before we continue with iteration "

# ╔═╡ b661b098-1d33-4470-8082-f3972e738e91
md" A closed form solution to this problem can be found if we choose $u(x) = \ln(x)$, but this will not always be the case. This is a very special type of problem that allows an analytic solution. "

# ╔═╡ 018646cc-912f-4ac9-9158-ebead586d663
begin
	# constructing the analytic solution to the problem. Useful for comparison wiht numerical answer.  
	ab        = α * β
	c₁        = (log(1 - ab) + log(ab) * ab / (1 - ab)) / (1 - β)
	c₂        = α / (1 - ab)
	
	# optimal analytical values
	v_star(k) = c₁ .+ c₂ .* log.(k)  
	k_star(k) = ab * k.^α   
	c_star(k) = (1-ab) * k.^α  
	ufun(x) = log.(x)
end

# ╔═╡ 828e7015-3b57-44fb-8cc4-0a7727245981
# Step 4: Start the iteration

# We will first define the Bellman operator before we can start the iteration. We will need this operator to perform the iteration. 

# Below is an example of what a docstring looks like for a function. 

"""
# Bellman Operator

# Inputs
`grid`: grid of values of state variable

`v0`: current guess of value function

`n` : number of grid points (default at `n₁` = 150 for now)

# Output
`v1`: next guess of value function

`pol`: corresponding policy function 
"""
function bellman_op(grid, v0, n, f; β = 0.96)
	
	# Initialise the different vectors
	v1  = zeros(n)     # next guess
    pol = zeros(Int,n)     # policy function
    w   = zeros(n)   # temporary vector 

	# we need to loop over the states (outer loop) and choices (inner loop)
    # loop over current states
    # current capital
    for (i,k) in enumerate(grid)

        # loop over all possible kprime choices
        for (iprime, kprime) in enumerate(grid)
            if f(k) - kprime < 0   #check for negative consumption
                w[iprime] = -Inf
            else
                w[iprime] = ufun(f(k) - kprime) + β * v0[iprime]
            end
        end
        # find maximal choice
        v1[i], pol[i] = findmax(w)     # stores value und policy (index of optimal choice) 
    end
    return (v1, pol)   # return both value and policy function
	
end

# ╔═╡ 4f944394-1175-4fc3-b3ee-b2185cab821b
b1 = bellman_op(K_grid, v₀, K_grid_size, k -> k^α)

# ╔═╡ 86869393-2a2e-4539-98b2-a061c90897c3
# VFI iterator
#
## input
# `n`: number of grid points
# output
# `v_next`: tuple with value and policy functions after `n` iterations.
function VFI(op)
    v_init = zeros(n)     # initial guess
    for iter in 1:N_iter
        v_next = op(kgrid,v_init)  # returns a tuple: (v1,pol)
        # check convergence
        if maximum(abs,v_init.-v_next[1]) < tol
            verrors = maximum(abs,v_next[1].-v_star(kgrid))
            perrors = maximum(abs,v_next[2].-k_star(kgrid))
            return (v = v_next[1], p =v_next[2], errv = verrors, errp = perrors, iter = iter)
        elseif iter==N_iter
            @warn "No solution found after $iter iterations"
            return (v = v_next[1], p =v_next[2], errv = verrors, errp = perrors, iter = iter)
        end
        v_init = v_next[1]  # update guess 
    end
end

# ╔═╡ da73f82c-8fe4-4aaf-8b0f-7e20cf27dfc6
function VFI_converge(op::Function,steps)
    v_init = zeros(n)     # initial guess
	pl = plot(kgrid, v_star(kgrid), color = :red, label = "true")
    for iter in 1:steps
		plot!(pl, kgrid, v_init, label = "", color = :grey)
        v_next = op(kgrid,v_init)  # returns a tuple: (v1,pol)
        # check convergence
        if maximum(abs,v_init.-v_next[1]) < tol
            verrors = maximum(abs,v_next[1].-v_star(kgrid))
            perrors = maximum(abs,v_next[2].-k_star(kgrid))
            return (v = v_next[1], p =v_next[2], errv = verrors, errp = perrors, iter = iter)
        elseif iter==N_iter
            @warn "No solution found after $iter iterations"
            return (v = v_next[1], p =v_next[2], errv = verrors, errp = perrors, iter = iter)
        end
        v_init = v_next[1]  # update guess 
    end
	pl
end

# ╔═╡ 4bc01bd5-17ef-42a7-b829-5bd369b938f3
md" #### Fitted value function iteration "

# ╔═╡ 2fba5422-a46d-4606-8c21-c0df54f3544b
md" The final topic for this session is fitted value function iteration. This method is quite similar to the value function iteration from before, but uses interpolation techniques. The general process looks like the following 

The process looks like this:

1. Begin with array of values $\{ V_0, \ldots, V_I \}$  representing the values of initial function $V$ on the grid points $\{ K_0, \ldots, K_I \}$.  
2. Build a function $\hat{V}$ on the state space $\mathbb{R}_+$ by linear interpolation, based on these data points.  
3. Obtain and record the value $T \hat{V}(K_i)$ on each grid point $K_i$  by repeatedly solving the maximization problem.  
4. Unless some stopping condition is satisfied, set $\{ V_0, \ldots, V_I \} = \{ T \hat V(K_0), \ldots, T \hat V(K_I) \}$ and go to step $2$.     "

# ╔═╡ 507514d4-d369-4950-9198-428561aecfd9
md" ##### Quick detour on interpolation "

# ╔═╡ ab73acc7-f51f-48cc-8b91-48d6d3a41c2b
md" We have covered the topic of function approximation via linear interapolation before, but let us have a quick recap before we continue. Suppose that we have a function $h(x) = 2\cos(6x) + \sin(14x) + 2.5$ and we are trying to approximate this function with piecewise linear interpolation, we could use the `Interpolations.jl` module. Consider the following:"

# ╔═╡ 1229edd5-45cf-41f1-ad66-72af83a0a2e3
begin
	CS1 = 0.2
	CS2 = 0.1
	CS3 = 0.05
end

# ╔═╡ 71fb6ee5-bd3f-4465-b2c6-e8e11c0865cb
@bind C_name Select(["CS1" => "5 Points", "CS2" => "10 Points", "CS3" => "20 Points"])

# ╔═╡ 1a3e85a7-a292-487e-804c-b408cd508498
c_size = Dict("CS1" => CS1, "CS2" => CS2, "CS3" => CS3)[C_name]

# ╔═╡ 597f6e04-3580-4214-b7f2-77db8286d864
begin
	h(x) = 2 .* cos.(6x) .+ sin.(14x) .+ 2.5
	c_grid = 0:c_size:1
	h_grid = range(0,  1, length = 150)
	
	Ah = LinearInterpolation(c_grid, h(c_grid))
	
	plt = plot(xlim = (0,1), ylim = (0,6))
	plot!(plt, h, h_grid, color = :black, lw = 2, alpha = 0.7, label = "true function")
	plot!(plt, h_grid, Ah.(h_grid), color = :green, lw = 2, alpha = 0.8,
	      label = "linear approximation")
	plot!(plt, h, c_grid, seriestype = :sticks, linestyle = :dash, linewidth = 2, alpha = 0.5,
	      label = "")
	plot!(plt, legend = :top)
end

# ╔═╡ fcca074c-8f19-42b3-b60e-a75cdad1311f
md" In this example we only discretize the state space, not the the action (control) space. We are providing an approximation to the control space via our interpolation. We mentioned this before, in the algorithm above. We are building a function $\hat{V}$ on $\mathbb{R}_{+}$. The values for $K^{\prime}$ are potentially then **off the grid** that we construct for the state space. "

# ╔═╡ 785151e1-f0b1-4788-833a-2db80120c876
num_grid = @bind m Slider(0.1:0.1:2,show_value = true, default = 0.2)

# ╔═╡ 4eea7313-37fe-48ea-b30a-94f4ccedd71d
scatter(K_grid, ones(1), xlims = (0, m), ylims=(0.9,1.1),label = "State space")

# ╔═╡ cbb90e8e-9656-4dfb-8103-6459fa46595e
md" For our Bellman operator we will use the `maximize()` function from the `Optim.jl` package to perform our maximisation. We define a objective function for each $k_{i}$ and then optimise with respect to objective, while defining a lower and upper bound. "

# ╔═╡ 2cef4d1b-39bd-4c24-8a3b-69c6da406e8f
md" The QuantEcon example utilises this exact process, but uses $y_t$ as state variable instead of $K_t$. This delivers exactly the same result, so we look at their code implementation below. In general, the notes and code at QuantEcon are much better thought out than I can manage, so I utilise their code whenever possible."

# ╔═╡ 91778068-9af3-4f28-bf76-aa19c9df3e45
function bellman_fitted_y(w, grid, β, u, f)
    w_func = LinearInterpolation(grid, w)
	
    # objective for each grid point
    objectives = (c -> u(c) + β * (w_func.(f(y - c))) for y in grid)
	
    results = maximize.(objectives, 1e-10, grid) # solver result for each grid point
  
	Tw = Optim.maximum.(results)
    return Tw
end

# ╔═╡ d111bb6f-f6ce-44bd-a40d-0d442d6484d0
function FVFI_y(w, grid, iter)
	
	lb = "initial condition"
	plt = plot(grid, w, color = :black, linewidth = 2, alpha = 0.8, label = lb)
	
	for i in 1:iter
		w = bellman_fitted_y(w, grid, β, log, k -> k^0.65)
		plot!(K_grid, w, color = RGBA(i/iter, 0.2, 1 - i/iter, 0.8), linewidth = 2, alpha = 0.7,
	          label = "")
	end
	plot!(plt, legend = :bottomright)
end

# ╔═╡ cac96345-ee00-40ef-b8ae-3b7fdbfea9be
w_new = 0.5 * log.(K_grid); 

# ╔═╡ 128ec2e2-3b38-49cd-b507-c1f5fb051491
iterations = @bind iter Slider(5:200, show_value = true, default = 5)

# ╔═╡ 07e571b5-9e5f-4690-89b6-4fc166c94a1b
FVFI_y(w_new, K_grid, iter)

# ╔═╡ d31741a8-5184-45d8-8a62-0db938c9a372
md" In the next session we will introduce a stochastic component to the problem, which means that there will be uncertainty in the model. In particular we will apply this to the cake eating and optimal growth models.  "

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ApproxFun = "28f2ccd6-bb30-5033-b560-165f7b14dc2f"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
Images = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
QuantEcon = "fcd29c91-0bd7-5a09-975d-7ac3f643a60c"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
ApproxFun = "~0.12.5"
Colors = "~0.12.8"
Images = "~0.24.1"
Interpolations = "~0.13.3"
LaTeXStrings = "~1.2.1"
NLsolve = "~4.5.1"
Optim = "~1.3.0"
Parameters = "~0.12.2"
Plots = "~1.18.2"
PlutoUI = "~0.7.9"
QuantEcon = "~0.16.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ApproxFun]]
deps = ["AbstractFFTs", "ApproxFunBase", "ApproxFunFourier", "ApproxFunOrthogonalPolynomials", "ApproxFunSingularities", "Calculus", "DomainSets", "DualNumbers", "FFTW", "FastTransforms", "LinearAlgebra", "RecipesBase", "Reexport", "SpecialFunctions"]
git-tree-sha1 = "411ca0706d880fb1fe173a9925477a3fe26c564f"
uuid = "28f2ccd6-bb30-5033-b560-165f7b14dc2f"
version = "0.12.5"

[[ApproxFunBase]]
deps = ["AbstractFFTs", "BandedMatrices", "BlockArrays", "BlockBandedMatrices", "Calculus", "DSP", "DomainSets", "DualNumbers", "FFTW", "FastGaussQuadrature", "FillArrays", "InfiniteArrays", "IntervalSets", "LazyArrays", "LinearAlgebra", "LowRankApprox", "SparseArrays", "SpecialFunctions", "StaticArrays", "Statistics", "Test", "ToeplitzMatrices"]
git-tree-sha1 = "bab00664565c82f38f562372951ca9b131072b5b"
uuid = "fbd15aa5-315a-5a7d-a8a4-24992e37be05"
version = "0.3.14"

[[ApproxFunFourier]]
deps = ["AbstractFFTs", "ApproxFunBase", "DomainSets", "FFTW", "FastTransforms", "InfiniteArrays", "IntervalSets", "LinearAlgebra", "Reexport"]
git-tree-sha1 = "50f8608bc46db8918ac80345bd033570b2e69bde"
uuid = "59844689-9c9d-51bf-9583-5b794ec66d30"
version = "0.2.9"

[[ApproxFunOrthogonalPolynomials]]
deps = ["AbstractFFTs", "ApproxFunBase", "BandedMatrices", "BlockArrays", "BlockBandedMatrices", "DomainSets", "FFTW", "FastGaussQuadrature", "FastTransforms", "FillArrays", "IntervalSets", "LinearAlgebra", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "a7329c5310e4c70e93516b96756c6cf5beafe905"
uuid = "b70543e2-c0d9-56b8-a290-0d4d6d4de211"
version = "0.3.10"

[[ApproxFunSingularities]]
deps = ["ApproxFunBase", "ApproxFunOrthogonalPolynomials", "DomainSets", "IntervalSets", "LinearAlgebra", "Reexport", "Statistics"]
git-tree-sha1 = "0d6851111d0cc741e845546bfc4a267a0f3a6433"
uuid = "f8fcb915-6b99-5be2-b79a-d6dbef8e6e7e"
version = "0.2.2"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[ArrayInterface]]
deps = ["IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "a71d224f61475b93c9e196e83c17c6ac4dedacfa"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.18"

[[ArrayLayouts]]
deps = ["Compat", "FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "8f6af27c33b766f19fa6cfe46e629775cda81f88"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.4.11"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "a4d07a1c313392a77042855df46c5f534076fab9"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.0"

[[AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "f31f50712cbdf40ee8287f0443b57503e34122ef"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.3"

[[BandedMatrices]]
deps = ["ArrayLayouts", "Compat", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "8c83ee44dc9835263760ad4e77ed4eed4b3490c1"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.15.25"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Statistics", "UUIDs"]
git-tree-sha1 = "c31ebabde28d102b602bada60ce8922c266d205b"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.1.1"

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
git-tree-sha1 = "c3598e525718abcc440f69cc6d5f60dda0a1b61e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.6+5"

[[CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "e2f47f6d8337369411569fd45ae5753ca10394c6"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.0+6"

[[Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f53ca8d41e4753c41cdafa6ec5f7ce914b34be54"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "0.10.13"

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
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random", "StaticArrays"]
git-tree-sha1 = "ed268efe58512df8c7e224d2e170afd76dd6a417"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.13.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "42a9b08d3f2f951c9b283ea427d96ed9f1f30343"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.5"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dc7dedc2c2aa9faf59a55c622760a25cbefbe941"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.31.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "6d1c23e740a586955645500bbec662476204a52c"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.1"

[[CustomUnitRanges]]
git-tree-sha1 = "537c988076d001469093945f3bd0b300b8d3a7f3"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.1"

[[DSP]]
deps = ["FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "2a63cb5fc0e8c1f0f139475ef94228c7441dc7d0"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.6.10"

[[DataAPI]]
git-tree-sha1 = "ee400abb2298bd13bfc3df1c412ed228061a2385"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.7.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4437b64df1e0adccc3e5d1adbc3ac741095e4677"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.9"

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

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "214c3fcac57755cfda163d91c58893a8723f93e9"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.0.2"

[[Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "abe4ad222b26af3337262b8afb28fab8d215e9f8"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.3"

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
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[DomainSets]]
deps = ["IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics", "Test"]
git-tree-sha1 = "372f681128afe9e728ff00dca2df6b0795e51fbe"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.4.2"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "fe385ec95ac5533650fb9b1ba7869e9bc28cdd0a"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.5"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "92d8f9f208637e8d2d28c664051a00569c01493d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.1.5+1"

[[EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "8041575f021cba5a099a456b4163c9a08b566a02"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.1.0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "LibVPX_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3cc57ad0a213808473eafef4845a74766242e05f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.3.1+4"

[[FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "70a0cfd9b1c86b0209e38fbfe6d8231fd606eeaf"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.1"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "f985af3b9f4e278b1d24434cbb546d6092fca661"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.3"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3676abafff7e4ff07bbd2c42b3d8201f31653dcc"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.9+8"

[[FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions"]
git-tree-sha1 = "6ea5f7b4aecce0e3a14ca1da03f62f86148c8fa3"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "0.4.5"

[[FastTransforms]]
deps = ["AbstractFFTs", "ArrayLayouts", "BinaryProvider", "DSP", "FFTW", "FastGaussQuadrature", "FastTransforms_jll", "FillArrays", "Libdl", "LinearAlgebra", "Reexport", "SpecialFunctions", "Test", "ToeplitzMatrices"]
git-tree-sha1 = "7c10535fb7b59bac11667749f2c1229ea755cad8"
uuid = "057dd010-8810-581a-b7be-e3fc3b93f78c"
version = "0.12.4"

[[FastTransforms_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "FFTW_jll", "JLLWrappers", "Libdl", "MPFR_jll", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "176f3f679f8921b3dc2ba127da2f9caf3f6a26eb"
uuid = "34b6f7d7-08f9-5794-9e10-3819e4c7e49a"
version = "0.5.1+0"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "256d8e6188f3f1ebfa1a5d17e072a0efafa8c5bf"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.10.1"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "502b3de6039d5b78c76118423858d981349f3823"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.9.7"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "35895cf184ceaab11fd778b4590144034a167a2f"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.1+14"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "NaNMath", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "e2af66012e08966366a43251e1fd421522908be6"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.18"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "cbd58c9deb1d304f5a245a0b7eb841a2560cfec6"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.1+5"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "dba1e8614e98949abfa60480b13653813d8f0157"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+0"

[[GMP_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "781609d7-10c4-51f6-84f2-b8444358ff6d"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "b83e3125048a9c3158cbb7ca423790c7b1b57bea"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.57.5"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "eaf96e05a880f3db5ded5a5a8a7817ecba3c7392"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.58.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "4136b8a5668341e58398bb472754bff4ba0456ff"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.3.12"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "47ce50b742921377301e15005c96e979574e130b"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.1+0"

[[Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "2c1cf4df419938ece72de17f368a021ee162762e"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "c6a1fff2fd4b1da29d3dccaffb1e1001244d844e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.12"

[[IdentityRanges]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be8fcd695c4da16a1d6d0cd213cb88090a150e3b"
uuid = "bbac6d45-d8f3-5730-bfe4-7a449cd117ca"
version = "0.3.1"

[[IfElse]]
git-tree-sha1 = "28e837ff3e7a6c3cdb252ce49fb412c8eb3caeef"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.0"

[[ImageAxes]]
deps = ["AxisArrays", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "794ad1d922c432082bc1aaa9fa8ffbd1fe74e621"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.9"

[[ImageContrastAdjustment]]
deps = ["ColorVectorSpace", "ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "2e6084db6cccab11fe0bc3e4130bd3d117092ed9"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.7"

[[ImageCore]]
deps = ["AbstractFFTs", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "db645f20b59f060d8cfae696bc9538d13fd86416"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.8.22"

[[ImageDistances]]
deps = ["ColorVectorSpace", "Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "6378c34a3c3a216235210d19b9f495ecfff2f85f"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.13"

[[ImageFiltering]]
deps = ["CatIndices", "ColorVectorSpace", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageCore", "LinearAlgebra", "OffsetArrays", "Requires", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "bf96839133212d3eff4a1c3a80c57abc7cfbf0ce"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.6.21"

[[ImageIO]]
deps = ["FileIO", "Netpbm", "PNGFiles", "TiffImages", "UUIDs"]
git-tree-sha1 = "d067570b4d4870a942b19d9ceacaea4fb39b69a1"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.5.6"

[[ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils", "Libdl", "Pkg", "Random"]
git-tree-sha1 = "5bc1cb62e0c5f1005868358db0692c994c3a13c6"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.2.1"

[[ImageMagick_jll]]
deps = ["JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "1c0a2295cca535fabaf2029062912591e9b61987"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.10-12+3"

[[ImageMetadata]]
deps = ["AxisArrays", "ColorVectorSpace", "ImageAxes", "ImageCore", "IndirectArrays"]
git-tree-sha1 = "ae76038347dc4edcdb06b541595268fca65b6a42"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.5"

[[ImageMorphology]]
deps = ["ColorVectorSpace", "ImageCore", "LinearAlgebra", "TiledIteration"]
git-tree-sha1 = "68e7cbcd7dfaa3c2f74b0a8ab3066f5de8f2b71d"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.2.11"

[[ImageQualityIndexes]]
deps = ["ColorVectorSpace", "ImageCore", "ImageDistances", "ImageFiltering", "OffsetArrays", "Statistics"]
git-tree-sha1 = "1198f85fa2481a3bb94bf937495ba1916f12b533"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.2.2"

[[ImageShow]]
deps = ["Base64", "FileIO", "ImageCore", "OffsetArrays", "Requires", "StackViews"]
git-tree-sha1 = "832abfd709fa436a562db47fd8e81377f72b01f9"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.1"

[[ImageTransformations]]
deps = ["AxisAlgorithms", "ColorVectorSpace", "CoordinateTransformations", "IdentityRanges", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "d966631de06f36c8cd4bec4bb2e8fa731db16ed9"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.8.12"

[[Images]]
deps = ["AxisArrays", "Base64", "ColorVectorSpace", "FileIO", "Graphics", "ImageAxes", "ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageShow", "ImageTransformations", "IndirectArrays", "OffsetArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "8b714d5e11c91a0d945717430ec20f9251af4bd2"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.24.1"

[[IndirectArrays]]
git-tree-sha1 = "c2a145a145dc03a7620af1444e0264ef907bd44f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "0.5.1"

[[InfiniteArrays]]
deps = ["DSP", "FillArrays", "LazyArrays", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "aa50bca542c26505efb399290eb5f83e995517bc"
uuid = "4858937d-0d70-526a-a4dd-2d5cb5dd786c"
version = "0.8.2"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

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
git-tree-sha1 = "1470c80592cf1f0a35566ee5e93c5f8221ebc33a"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.3"

[[IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "81690084b6198a2e1da36fcfda16eeca9f9f24e4"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.1"

[[JSONSchema]]
deps = ["HTTP", "JSON", "ZipFile"]
git-tree-sha1 = "b84ab8139afde82c7c65ba2b792fe12e01dd7307"
uuid = "7d188eb4-7ad8-530c-ae41-71a32a6d4692"
version = "0.3.3"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

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

[[LibVPX_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "12ee7e23fa4d18361e7c2cde8f8337d4c3101bc7"
uuid = "dd192d2f-8180-539f-9fb4-cc70b1dcf69a"
version = "1.10.0+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "761a393aeccd6aa92ec3515e428c26bf99575b3b"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+0"

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
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

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
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["DocStringExtensions", "LinearAlgebra"]
git-tree-sha1 = "7bd5f6565d80b6bf753738d2bc40a5dfea072070"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.2.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LowRankApprox]]
deps = ["FFTW", "FillArrays", "LinearAlgebra", "Nullables", "Random", "SparseArrays", "Test"]
git-tree-sha1 = "35e31d7e505492cb70b44eacf94e89adf6ef79f6"
uuid = "898213cb-b102-5a47-900c-97e73b919f73"
version = "0.4.3"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "c253236b0ed414624b083e6b72bfe891fbd2c7af"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+1"

[[MPFR_jll]]
deps = ["Artifacts", "GMP_jll", "Libdl"]
uuid = "3a97d323-0669-5f0c-9066-3539efd106a3"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6a8a2a625ab0dea913aba95c11370589e0239ff0"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.6"

[[MappedArrays]]
git-tree-sha1 = "18d3584eebc861e311a552cbb67723af8edff5de"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.0"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "JSON", "JSONSchema", "LinearAlgebra", "MutableArithmetics", "OrderedCollections", "SparseArrays", "Test", "Unicode"]
git-tree-sha1 = "575644e3c05b258250bb599e57cf73bbf1062901"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "0.9.22"

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
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "4ea90bd5d3985ae1f9a908bd4500ae88921c5ce7"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.0"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["ExprTools"]
git-tree-sha1 = "916b850daad0d46b8c71f65f719c49957e9513ed"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.1"

[[MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3927848ccebcc165952dc0d9ac9aa274a87bfe01"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.20"

[[NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50608f411a1e178e0129eab4110bd56efd08816f"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.0"

[[NLopt]]
deps = ["MathOptInterface", "MathProgBase", "NLopt_jll"]
git-tree-sha1 = "d80cb3327d1aeef0f59eacf225e000f86e4eee0a"
uuid = "76087f3c-5699-56af-9a33-bf431cd00edd"
version = "0.6.3"

[[NLopt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "2b597c46900f5f811bec31f0dcc88b45744a2a09"
uuid = "079eb43e-fd8e-5478-9966-2cf3e3edb778"
version = "2.7.0+0"

[[NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[Netpbm]]
deps = ["ColorVectorSpace", "FileIO", "ImageCore"]
git-tree-sha1 = "09589171688f0039f13ebe0fdcc7288f50228b52"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.1"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Nullables]]
git-tree-sha1 = "8f87854cc8f3685a60689d8edecaa29d2251979b"
uuid = "4d1e1d77-625e-5b40-9113-a560ec7a8ecd"
version = "1.0.0"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "4f825c6da64aebaa22cc058ecfceed1ab9af1c7e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.3"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Optim]]
deps = ["Compat", "FillArrays", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "d34366a3abc25c41f88820762ef7dfdfe9306711"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.3.0"

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

[[PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "520e28d4026d16dcf7b8c8140a3041f0e20a9ca8"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.7"

[[PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fa5e78929aebc3f6b56e1a88cf505bb00a354c4"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.8"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "2276ac65f1e236e0a6ea70baff3f62ad4c625345"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.2"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "c8abc88faa3f7a3950832ac5d6e690881590d6dc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.0"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "501c20a63a34ac1d015d5304da0e645f42d91c9f"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.11"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "f32cd6fcd2909c2d1cdd47ce55e1394b04a66fe2"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.18.2"

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
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Primes]]
git-tree-sha1 = "afccf037da52fa596223e5a0e331ff752e0e845c"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "afadeba63d90ff223a6a48d2009434ecee2ec9e8"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.1"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "12fbe86da16df6679be7521dfb39fbc861e1dc7b"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.1"

[[QuantEcon]]
deps = ["DSP", "DataStructures", "Distributions", "FFTW", "LightGraphs", "LinearAlgebra", "Markdown", "NLopt", "Optim", "Pkg", "Primes", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "Test"]
git-tree-sha1 = "4e2dc3044303aa2cbf6e321cb9af3982f6774e6a"
uuid = "fcd29c91-0bd7-5a09-975d-7ac3f643a60c"
version = "0.16.2"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[Ratios]]
git-tree-sha1 = "37d210f612d70f3f7d57d488cb3b6eff56ad4e41"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.0"

[[RecipesBase]]
git-tree-sha1 = "b3fb709f3c97bfc6e948be68beeecb55a0b340ae"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "2a7a2469ed5d94a98dea0e85c46fa653d76be0cd"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.3.4"

[[Reexport]]
git-tree-sha1 = "5f6c21241f0f655da3952fd60aa18477cf96c220"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.1.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

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

[[Rotations]]
deps = ["LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "2ed8d8a16d703f900168822d83699b8c3c1a5cd8"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.0.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

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

[[StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "62701892d172a2fa41a1f829f66d2b0db94a9a63"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.3.0"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "da4cf579416c81994afd6322365d00916c79b8ae"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "0.12.5"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2f6792d523d7448bbe2fec99eca9218f06cc746d"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.8"

[[StatsFuns]]
deps = ["LogExpFunctions", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "30cd8c360c54081f806b1ee14d2eecbef3c04c49"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.8"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "Tables"]
git-tree-sha1 = "44b3afd37b17422a62aea25f04c1f7e09ce6b07f"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.5.1"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "8ed4a3ea724dac32670b062be3ef1c1de6773ae8"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.4.4"

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

[[TiffImages]]
deps = ["ColorTypes", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "OffsetArrays", "OrderedCollections", "PkgVersion", "ProgressMeter"]
git-tree-sha1 = "03fb246ac6e6b7cb7abac3b3302447d55b43270e"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.4.1"

[[TiledIteration]]
deps = ["OffsetArrays"]
git-tree-sha1 = "52c5f816857bfb3291c7d25420b1f4aca0a74d18"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.3.0"

[[TimeZones]]
deps = ["Dates", "Future", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "81753f400872e5074768c9a77d4c44e70d409ef0"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.5.6"

[[ToeplitzMatrices]]
deps = ["AbstractFFTs", "FFTW", "LinearAlgebra", "StatsBase"]
git-tree-sha1 = "39731ea29652305b1eeb7ab0d1085533db307191"
uuid = "c751599d-da0a-543b-9d20-d0a503d91d24"
version = "0.6.3"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "7c53c35547de1c5b9d46a4797cf6d8253807108c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.5"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "59e2ad8fd1591ea019a5259bd012d7aee15f995c"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.3"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

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

[[ZipFile]]
deps = ["Libdl", "Printf", "Zlib_jll"]
git-tree-sha1 = "c3a5637e27e914a7a445b8d0ad063d701931e9f7"
uuid = "a5390f91-8eb1-5f08-bee0-b1d1ffed6cea"
version = "0.9.3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "acc685bcf777b2202a904cdcb49ad34c2fa1880c"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.14.0+4"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7a5780a0d9c6864184b3a2eeeb833a0c871f00ab"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "0.1.6+4"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d713c1ce4deac133e3334ee12f4adff07f81778f"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2020.7.14+2"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "487da2f8f2f0c8ee0e83f39d13037d6bbf0a45ab"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.0.0+3"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─7fffd413-6ebe-49d2-b426-1428aee5ae26
# ╟─c428b2b9-a4f4-4ddc-be20-d25cb14c9cf7
# ╟─56121f90-e33f-11eb-1211-8578ae4eb05d
# ╟─058c997b-390b-4b8e-b203-20879e0673a4
# ╟─f6852b65-d493-4e2c-a8b6-517993c93ac8
# ╟─3e79f184-a3e8-44a0-9609-19827522d86a
# ╟─340698bc-e27a-442b-a6fe-4b7c7ed3bd16
# ╟─08bf10d7-0221-40e8-99ac-4087ad459add
# ╟─6071a006-d300-487b-835e-bc4b8487d458
# ╟─621779a2-fcbd-4e55-8b36-447733b92ea5
# ╟─eea67d06-ac30-4e49-b307-b6b4691f1d76
# ╟─f6c5bc76-b0d8-46e0-9994-b6163eee9704
# ╟─fd77c87a-ac60-4507-a251-19cb4b82088a
# ╟─0c7a94ba-dd82-4015-9c39-adaec2e424f5
# ╟─d2de514f-ae9d-4d89-ba5e-e5a2d4c55081
# ╟─1fa1bdc2-a076-42be-9540-1d20ebcffeb5
# ╟─e5e5bf22-a1d3-4bc3-8968-aa1d0545c2c8
# ╟─b844d7eb-b152-4c37-9922-a0223f3b5589
# ╟─ba5f4d25-2ed0-4d5b-8f84-25f10dc0ad82
# ╟─e5ef57d1-e00b-4005-913a-7d6579c9c59c
# ╟─398a316b-fc53-444a-aea7-169231c95f4b
# ╟─b0ca2d8d-dfa9-4fed-b34d-90167e950153
# ╟─db4c8010-a33f-4e78-bea4-b3168c5ce72b
# ╟─246e0451-111d-423c-a76f-eb4339b2d285
# ╟─728e0ddc-f482-4d9a-b6f2-41f7cf6b8619
# ╟─3f18f963-ad85-4c75-ab8c-64ac3cbf312c
# ╟─fc6640d2-f2de-4d09-8787-39954cc835b4
# ╟─a897d516-60e9-4b96-9f1f-f7dea164f392
# ╟─f17bd25f-a459-44b6-862e-5d8d7b5062c9
# ╟─be62b6ff-a947-47ab-8905-6ac631abf8e6
# ╟─cd9f3fd5-307a-4e1d-b92f-1bf7d441d483
# ╠═d9a659a2-d0dd-45bb-bb5b-fb85f0208819
# ╠═a25b4e49-46c0-4c7d-b80e-e108614ec288
# ╠═11dc1b7f-866e-42e6-b4ef-0b9d3b34366f
# ╟─f2075369-02b7-4262-a543-91da85169e65
# ╟─8b8c03cb-b8c5-45a8-9d45-a6730f0cf420
# ╠═ed53fb3c-9fe6-47d5-b465-3af3369f9177
# ╠═d609467a-6898-430f-9fc7-70cecf9cb3d0
# ╠═be24ede3-87c7-402c-adbc-2f003e5099a0
# ╠═8e14ac14-af81-4bd6-bdac-f850bb9020e9
# ╠═9300d7d9-9563-4f85-b835-54a8db884317
# ╟─53bdea85-8d59-4669-85ca-04119131064a
# ╟─eff5cdcf-3a07-4d08-9152-39c0b71eac3b
# ╠═b6cbcf47-669d-4ae3-aa79-b83d59e04d9d
# ╟─437462d4-2dd1-4ef9-aa7d-b0cb64bef455
# ╠═18b69513-700f-4300-b578-4742112a23ca
# ╟─53c21b1b-088f-4604-9246-fc6796c3e256
# ╟─fd4d1601-5a71-4ce1-a3bd-5d34776502f9
# ╟─c8533761-ae1d-47a2-9860-ad88ecca2d34
# ╟─547def39-8a65-48bc-bb30-d28d016cd6fa
# ╟─7b63094f-ddde-4dff-911a-58d8fc7226a2
# ╟─c5ea0324-8e8c-49d1-88f6-3abc4278aefa
# ╠═4716e611-3342-4c3e-afc2-5afd8de2fc15
# ╟─aa534428-e995-4885-bfeb-135ad92129a4
# ╟─b40713e0-54db-4b49-8d35-a5b68f71ac89
# ╟─7854bba6-88c3-4f19-b7b1-1ae8487ac9bc
# ╟─b71abd7b-191f-4044-91a9-d0f4e4883f88
# ╟─6f2eac77-4ad3-4da7-95a7-c958ccf270cd
# ╟─d81d8f8a-ac8a-4d01-89c8-bb62960a69bc
# ╟─7529e1a7-c41c-49ec-8885-6c7c14093700
# ╟─b22748fb-9acb-4105-9fbe-654daf34dbf4
# ╠═86f2d8ca-ea32-4558-8133-bce4789ad105
# ╟─b1e0e8e1-858f-4f12-a272-7197a117af25
# ╟─23bdbe50-8544-4bb9-9e91-ba7397ca4db2
# ╟─8ec82051-79b0-4cd7-8585-163ffde2b290
# ╟─616d4193-6ba7-450e-ad3e-d6b07cb90d2e
# ╟─bfeff134-00d8-4ca5-9902-6b8dde5facd3
# ╟─f4654f14-6bef-452f-a9af-ec62613ad20b
# ╟─5a3df38d-2973-421e-8da0-a190653af20f
# ╟─e8e3eb80-0c0f-4a08-b963-eb554cac8f0b
# ╟─8f7abbcf-61c5-4c4f-92f4-41c02a6aaa41
# ╟─51c0c398-81d0-4e93-9907-34f3f39be2e4
# ╟─b2f96156-bd23-4767-8ec8-38215a8f4b63
# ╠═60912ad4-c766-4b70-b843-baed710e29ce
# ╠═adc503f5-6f65-4f6a-ba10-965f78f52e9e
# ╠═cd71f975-5058-46a4-80d8-0d410f19cbf1
# ╟─a8cbe356-218d-47f1-b62a-8fea4f6c1835
# ╠═828e7015-3b57-44fb-8cc4-0a7727245981
# ╠═4f944394-1175-4fc3-b3ee-b2185cab821b
# ╠═86869393-2a2e-4539-98b2-a061c90897c3
# ╠═da73f82c-8fe4-4aaf-8b0f-7e20cf27dfc6
# ╟─b661b098-1d33-4470-8082-f3972e738e91
# ╠═018646cc-912f-4ac9-9158-ebead586d663
# ╟─4bc01bd5-17ef-42a7-b829-5bd369b938f3
# ╟─2fba5422-a46d-4606-8c21-c0df54f3544b
# ╟─507514d4-d369-4950-9198-428561aecfd9
# ╟─ab73acc7-f51f-48cc-8b91-48d6d3a41c2b
# ╟─1229edd5-45cf-41f1-ad66-72af83a0a2e3
# ╟─71fb6ee5-bd3f-4465-b2c6-e8e11c0865cb
# ╟─1a3e85a7-a292-487e-804c-b408cd508498
# ╟─597f6e04-3580-4214-b7f2-77db8286d864
# ╟─fcca074c-8f19-42b3-b60e-a75cdad1311f
# ╟─785151e1-f0b1-4788-833a-2db80120c876
# ╠═4eea7313-37fe-48ea-b30a-94f4ccedd71d
# ╟─cbb90e8e-9656-4dfb-8103-6459fa46595e
# ╟─2cef4d1b-39bd-4c24-8a3b-69c6da406e8f
# ╠═91778068-9af3-4f28-bf76-aa19c9df3e45
# ╠═d111bb6f-f6ce-44bd-a40d-0d442d6484d0
# ╠═cac96345-ee00-40ef-b8ae-3b7fdbfea9be
# ╟─128ec2e2-3b38-49cd-b507-c1f5fb051491
# ╠═07e571b5-9e5f-4690-89b6-4fc166c94a1b
# ╟─d31741a8-5184-45d8-8a62-0db938c9a372
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
