### A Pluto.jl notebook ###
# v0.12.12

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

# ╔═╡ 5b44b6ec-7863-11eb-1ed4-3d0c9eadd065
using Markdown

# ╔═╡ 66ca1514-7863-11eb-2fd9-c5a2e7d1a74d
using InteractiveUtils

# ╔═╡ 6f5370c6-7b78-11eb-2076-45dc8e4908ae
using GraphRecipes, Plots

# ╔═╡ d482635e-7ab4-11eb-0d53-398ca84dab54
using DecFP

# ╔═╡ 8598d0ca-7863-11eb-1549-21bc81a5cb1f
md"# Fundamentals of numerical methods

The goal of this session is to provide a brief overview of the elementary ideas used in numerical methods and approximation procedures. The discussion in this section mostly focusses on sources of inexactness in computational science, which are referred to as approximations. Errors generally arise from two basic aspects of numerical operations, rounding and truncation. We will talk about truncation errors at another time, for now the focus is on rounding errors. In order to appreciate rounding errors, we need to better understand floating point numbers and computer arithmetic. "

# ╔═╡ e6a2f610-7ac6-11eb-2920-0ba6b5af0ee6
md" ### Floating point numbers "

# ╔═╡ f1f991e0-7864-11eb-1b11-abf2aa055858
md"In this section we introduce the idea of floating point numbers. It is not crucial for the rest of the work, but it is useful to understand. Here is a nice video that motivates the usage of floating point numbers. I like the informal way in which it is discussed in the video."

# ╔═╡ e25b1a6a-7b75-11eb-06de-bdfec783caa1
html"""

<iframe width="892" height="502" src="https://www.youtube.com/embed/PZRI1IfStY0" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

"""


# ╔═╡ e6015b5c-7b75-11eb-2671-b19f7e366107
md"You might have heard the term `binary` before. This refers to a base 2 number system. Each digit can take on the value of zero or one. `bit` stands for `binary` + `digit`. A `byte` is equal to 8 `bits`. From there we have that a `kilobyte` is $$10^3$$ `bytes` and a `megabyte` is $$10^6$$ `bytes`. This pattern continues for gigabytes, terabytes, petabytes, exabytes, zetabytes and so forth."

# ╔═╡ 85552d00-7865-11eb-1a94-3d5d52ed1772
md" As an aside, we can determine the amount of memory used by an object with the `Base.summarysize` command in _Julia_. I am not aware a specific function that does this in _Matlab_, there is a `whos()` funtion that returns the size of all variables in the workspace. The _Julia_ equivalent of `whos()` is `varinfo()`."

# ╔═╡ 84df01a6-7861-11eb-2e3d-bd2e48f379f0
x = rand(10, 10)

# ╔═╡ 6db2ed2a-7862-11eb-0053-99ca2c0cc1e6
Base.summarysize(x)

# ╔═╡ 305a2bfe-7863-11eb-06c7-39c0544baece
md"A fixed point number system is a computer model for integers $$\mathbb{Z}$$, while the floating number system is a computer model for real numbers $$\mathbb{R}$$. We want to essentially use a finite set to replace real numbers. The fact that the real numbers are infinite along to dimensions (continuous and uncountable) is problematic in computing. One solution is to use the floating point numbers. "

# ╔═╡ 3d20faae-78c5-11eb-3250-0932a784a744
md"The floating point representation of a real number is $$\pm d_0.d_1d_2 \ldots d_p \times b^e$$. For the computer the base is $$b = 2$$ and the digits $$d_i$$ are zero or one.

In the floating number system the computer stores a `sign bit`, the `significand` and the actual `exponent bit` plus a bias.  

The first bit is a `sign bit`. Then there are $11$ `exponent bits`. Finally, we have $p = 52$ bits that form the `significand`. The $52$ bit `significand` precision gives $15$ significant decimal digit precision.  

The following image is retrieved from the Wipedia entry on double floating point format and helps visualise the floating point representation."

# ╔═╡ 7f430c70-7b4e-11eb-3f09-7bd60114d1de
md" ![floating point numbers](https://upload.wikimedia.org/wikipedia/commons/thumb/a/a9/IEEE_754_Double_Floating_Point_Format.svg/1024px-IEEE_754_Double_Floating_Point_Format.svg.png)"

# ╔═╡ ed6adbd2-7b4a-11eb-2ec3-b7f98e04c983
md"From this picture you can see which parts form the `sign`, `exponent` and `significand`.  "

# ╔═╡ d5191c50-7b55-11eb-20c9-896d13cd253e
md"In _Julia_ (and most other programming languages) `1` and `1.0` are different values, because they have different types. **Types** are incredibly important and we will talk more about them at a later stage. We can determine the type of any object by using the `typeof()` command."

# ╔═╡ ee60791a-7b55-11eb-0bed-2d54e6d2f18b
typeof(1)

# ╔═╡ f92a2226-7b55-11eb-2dd1-bbdb6d228fef
typeof(1.0)

# ╔═╡ a0592f30-7b78-11eb-2482-f5eaa3a66228
md" I have plotted the type hierarchy for floating point numbers below. Don't worry if you don't understand this, I will explain in the meeting what this means."

# ╔═╡ 65b2d784-7b78-11eb-063b-1b7c30edd211
begin
	#pyplot(size=(800, 600))
	gr(size=(1000, 1000))
	theme(:default)
	
	plot(AbstractFloat, method=:tree, fontsize=8, nodeshape=:ellipse)
end

# ╔═╡ 4e406a14-7b56-11eb-109c-6136655fc096
md"Double precision numbers are the dominant data type in scientific computing. `Float64` is the type for double precision numbers in _Julia_. As far as I know, _Matlab_ only uses the double precision floating point number. These numbers utilise $64$ bits of memory. We can use the `bitstring()` command to see the bitstring representation."

# ╔═╡ fd8a93ee-7b55-11eb-0574-5fe37395446f
bitstring(1.0)

# ╔═╡ fa7f40a4-7b56-11eb-196f-81682455b932
md" We can access the different components of the floating point values as follows:"

# ╔═╡ e4769b72-7b56-11eb-3a0f-6f19f7325ab1
g = 1.0; @show g, sign(g), exponent(g), significand(g)

# ╔═╡ c29d3700-7b55-11eb-16ce-41d7accc401a
md"In the example below we use a slider to give an idea of what happens as we change the value of the number. Shift the slider to see what happens with the bitstring representation."

# ╔═╡ 400f0d96-7b45-11eb-1aec-df466a986bed
@bind y Slider(-9.9:0.1:9.9, default=1)

# ╔═╡ 535617f8-7b49-11eb-320e-8f7517146820
y

# ╔═╡ a12171e6-7ac8-11eb-1961-057baa390013
typeof(y)

# ╔═╡ 1fddcac8-7b48-11eb-0ed3-51b5b644e675
bitstring(y)

# ╔═╡ 5f54346a-7866-11eb-3df4-658de2c2ce5f
md"### Floating point arithmetic"

# ╔═╡ 5cf34e4e-7867-11eb-18ff-01b71b0d5934
md"From the previous section we know that all numbers in a computer are represented in binary format (zeros and ones). This can lead to some peculiar arithmetic properties even when you have simple expressions. It is super important to know how computers work with their representation of *real numbers*, otherwise we will make many mistakes. Let us show a simple example to illustrate."

# ╔═╡ ac567efc-7867-11eb-07c3-21989feec15c
begin
	a = 0.1
	b = 0.1
	c = 0.1
	
	a + b + c == 0.3
end;

# ╔═╡ 31795f9a-7b7a-11eb-1786-89aad77aff4b
md" The reason that we have this issue is due to rounding. Rounding is necesarry when a number has more than $p$ `significand` bits. _Julia_ offers several rounding modes, with the default being `RoundNearest`. In our example here, the number `0.1` in the decimal system cannot be represented accurately as a floating point number." 

# ╔═╡ de5fc23a-7b7a-11eb-12b0-a513044b39a6
bitstring(0.1)

# ╔═╡ e7475a98-7b7a-11eb-00a2-63e52e403c16
md"Do you see the repeating pattern in the `significand`? Even thought `0.1` can easily be represented with only two base-10 digits, it requires an infinite number of binary digits, which is cut off."

# ╔═╡ 3310e738-7ac1-11eb-38ea-2f1cfc2f0090
md"Let us illustrate with another example that exact arithmetic and computer arithmetic don't always give the same answers. Consider the following:"

# ╔═╡ 8f365958-7ac1-11eb-145b-69bff762b7a8
begin 
	e = (1e-20 + 1) - 1
	f = 1e-20 + (1 - 1)
	e == f
end;

# ╔═╡ 9fb90230-7ac1-11eb-1c50-634ba7a9a7bb
md" The value we obtain for `e` is actually incorrect."

# ╔═╡ f3fce35a-7b77-11eb-0935-13399babee72
e;

# ╔═╡ f80d6c4e-7b77-11eb-31bc-d55c9905b0a6
f;

# ╔═╡ fbf35d0a-7b77-11eb-3b36-194f011f00ab
md" Even though the statements are mathematically the same (by the associative property of addition and subtraction), if we compare the values we see that `e` > `f`. In this case this is because adding numbers of different magnitudes does not always work like you would want. In the example, when we added $10^{-20}$ to $1$, it got rounded away. "

# ╔═╡ 0b370ade-7b5d-11eb-14ae-8fb84aa69edc
md"Another class of problems using floating point arithmentic is illustrated in the following example:"

# ╔═╡ 19fd463c-7b5d-11eb-342c-293f27a0f396
begin 
	h = 100000.2
	i = 100000.1
end;

# ╔═╡ 2f54327a-7b5d-11eb-0352-e713aba2460d
h - i;

# ╔═╡ 3a62156a-7b5d-11eb-1b66-a9d7d71d1051
md"Think about what you would expect the answer to the calculation to be. Did you get it right?"

# ╔═╡ 96d3fed8-7b7b-11eb-1434-07da91c27aca
md" What about the following example. Think about what the result should be and what the computer calculates. Did you get this one right? "

# ╔═╡ c5f3dde8-7b7b-11eb-3a14-bbe7f374c856
begin
	j = 2.0^30 # Very large number
	k = 2.0^-30 # Very small number
	j + k == j
end;

# ╔═╡ 2a044f08-7b7b-11eb-19fc-9fe29ad589bf
md"Both of the examples above showcase the idea of **catastrophic cancellation**. In the first instance the subtraction of two nearly equal numbers eleminates significant digits (since the numbers can't be precisely represented by the machine) and in the second case there is a potential loss of precision when subtracting close real numbers that are represented by floats. 

For the first example we see that the values are approximately equal, as shown below."

# ╔═╡ 8b29e5b4-7b7f-11eb-2679-bd248ed2d88b
isapprox(100000.2 - 100000.1, 0.1)

# ╔═╡ 0f9373f0-78c5-11eb-26ba-05bea8874141
md" These examples makes it somewhat difficult for us to trust calculations if we don't know what the pitfalls are. Look at the following investment calculation. Can we trust this answer?"

# ╔═╡ 3933c4ee-78c5-11eb-2692-7989afac7ebb
begin
	interest = 0.04 # Interest rate
	compounding = 365 * 24 # Compounded every hour for 365 days
	investment = 10e9 # Initial investment of R1 billion
	t = 100 # Time invested (100 years)
	daily = 1 + interest / compounding # Daily return
	sum = investment * (daily ^ (compounding * t)) # Simple investment calculation
end;

# ╔═╡ 6dda76f6-78c6-11eb-0fe0-a36a5663dfd7
md"We have done the calculation using floating point numbers, now we use the precise decimal representation for a comparison."

# ╔═╡ add65ba2-78c7-11eb-0eb6-771466c0c417
import Pkg; Pkg.add("DecFP") 

# ╔═╡ 844fa292-7b48-11eb-0f9c-c521a4125ce6
Pkg.add("PlutoUI"); using PlutoUI

# ╔═╡ 70ab6652-7ab4-11eb-1309-8f833821ef88
begin
	interest_1 = Dec64(interest)
	daily_1 = 1 + interest_1 / compounding
	investment_1 = Dec64(investment)
	sum_1 = investment_1*(daily_1 ^ (compounding * t)) # Using exact decimals
end;

# ╔═╡ 3cf74998-78c5-11eb-1191-a9f56dbadda6
md"Compare the results and see if there is a significant difference. "

# ╔═╡ 3cfda932-78c5-11eb-3ad9-512459ffebc4
diff = sum_1 - Dec64(sum);

# ╔═╡ 3d11b968-78c5-11eb-3962-f50b8e05803d
md" Do you expect there to be any difference between the answers?"

# ╔═╡ 3d198bac-78c5-11eb-35af-87b99956e4de
md"We have to realise that real numbers are only presented with a certain precision. Most of the time we don't need to worry, but we want to make sure our code is robust. In other words, we shouldn't be surprised by the outcome of our code. The **numerical stability** of our code is important! "

# ╔═╡ 9d74de0e-7ac6-11eb-10c2-e3ac9cf17111
md" One could also have **overflow** (**underflow**) problems. In these instances, this entails a real number that is too large (indistinguishable from zero)."

# ╔═╡ 7bd38294-7b7c-11eb-0329-99e46b366e1b
md" ## Rules for numerical computation "

# ╔═╡ 64484bfa-7b7c-11eb-1d40-296fcf167f53
md" Some basic rules that we want to follow when doing computation.

1. Add small numbers together before adding larger ones
2. Add numbers of like magnitude together. This is called pairing.
3. Avoid subtraction of two numbers that are nearly equal. (Subtractive cancellation)"

# ╔═╡ 0ff8fb0c-7b84-11eb-2a6a-0dedd44b3c7b
md" Note that subtractive cancellation is one of the most common mechanisms introducing dramatic growth of errors in floating point computation "

# ╔═╡ cc6069f2-7b7c-11eb-3cb9-afbd773708e2
md" You could have saved a lot of time and just looked at these rules, but the examples give an idea as to why you might want to implement the rules in the first place."

# ╔═╡ 18b49d34-7b80-11eb-3b12-3774bdf47ae6
md" We will encounter more types of errors in approximation in the future, especially when we talk about calculating derivatives. Then we will introduce ways in which we can quantify the error with different metrics. However, I think this is enough information for one session. Your next step is to look at some of the exercises below."

# ╔═╡ b2b0b842-7b57-11eb-1727-95c07b5b2de6
md" ## Exercises 

For the first few exercises let's check whether floating-point numbers obey certain algebraic rules. For 2-5, one couterexample will suffice. 

1. The associative rule for addition says `(x + y) + z == x + (y + z)`. Check this rule using `x = 0.1`, `y = 0.1` and `z = 1.0`. Explain what you find.
2. Do floating-point numbers obey the associative rule for multiplication: `(x * y) * z == x * (y * z)`?
3. Do floating-point numbers obey the distributive rule: `a * (x + y) == a * x + a * y`?
4. Is `0 * x == 0` true for all `x` (where `x` is a floating point number)?
5. Is `x / a == x * (1 / a)` always true?"

# ╔═╡ 7779ed4c-7b7d-11eb-034b-5fe2925df5df
md" 6. **[Hard]** In this problem we will use `for loops`, so be sure that you have an idea of how to use loops. We can expect floating point errors to accumulate randomly during computation, creating what is known as a **random walk**. On average we expect the errors to partially cancel out. Suppose you define a random sequence by $x_0 = 0$ and $x_n = x_{n-1} \pm 1$ with the signs chosen by tossing a fair coin for each $n$. et $\alpha_n$ and $\beta_n$ be the average value of $x_n$ and $|x_n|$ respectively over all such walks. Then a classic result of probability is that $\alpha_n = 0$ and $\lim_{n\rightarrow \infty} \frac{\pi \beta^{2}_{n}}{2n} = 1$. Perform a million random walks"

# ╔═╡ Cell order:
# ╟─8598d0ca-7863-11eb-1549-21bc81a5cb1f
# ╠═5b44b6ec-7863-11eb-1ed4-3d0c9eadd065
# ╠═66ca1514-7863-11eb-2fd9-c5a2e7d1a74d
# ╠═844fa292-7b48-11eb-0f9c-c521a4125ce6
# ╟─e6a2f610-7ac6-11eb-2920-0ba6b5af0ee6
# ╟─f1f991e0-7864-11eb-1b11-abf2aa055858
# ╟─e25b1a6a-7b75-11eb-06de-bdfec783caa1
# ╟─e6015b5c-7b75-11eb-2671-b19f7e366107
# ╟─85552d00-7865-11eb-1a94-3d5d52ed1772
# ╠═84df01a6-7861-11eb-2e3d-bd2e48f379f0
# ╠═6db2ed2a-7862-11eb-0053-99ca2c0cc1e6
# ╟─305a2bfe-7863-11eb-06c7-39c0544baece
# ╟─3d20faae-78c5-11eb-3250-0932a784a744
# ╟─7f430c70-7b4e-11eb-3f09-7bd60114d1de
# ╟─ed6adbd2-7b4a-11eb-2ec3-b7f98e04c983
# ╟─d5191c50-7b55-11eb-20c9-896d13cd253e
# ╠═ee60791a-7b55-11eb-0bed-2d54e6d2f18b
# ╠═f92a2226-7b55-11eb-2dd1-bbdb6d228fef
# ╟─a0592f30-7b78-11eb-2482-f5eaa3a66228
# ╟─6f5370c6-7b78-11eb-2076-45dc8e4908ae
# ╟─65b2d784-7b78-11eb-063b-1b7c30edd211
# ╟─4e406a14-7b56-11eb-109c-6136655fc096
# ╠═fd8a93ee-7b55-11eb-0574-5fe37395446f
# ╟─fa7f40a4-7b56-11eb-196f-81682455b932
# ╠═e4769b72-7b56-11eb-3a0f-6f19f7325ab1
# ╟─c29d3700-7b55-11eb-16ce-41d7accc401a
# ╟─400f0d96-7b45-11eb-1aec-df466a986bed
# ╟─535617f8-7b49-11eb-320e-8f7517146820
# ╟─a12171e6-7ac8-11eb-1961-057baa390013
# ╟─1fddcac8-7b48-11eb-0ed3-51b5b644e675
# ╟─5f54346a-7866-11eb-3df4-658de2c2ce5f
# ╟─5cf34e4e-7867-11eb-18ff-01b71b0d5934
# ╠═ac567efc-7867-11eb-07c3-21989feec15c
# ╟─31795f9a-7b7a-11eb-1786-89aad77aff4b
# ╠═de5fc23a-7b7a-11eb-12b0-a513044b39a6
# ╟─e7475a98-7b7a-11eb-00a2-63e52e403c16
# ╟─3310e738-7ac1-11eb-38ea-2f1cfc2f0090
# ╠═8f365958-7ac1-11eb-145b-69bff762b7a8
# ╟─9fb90230-7ac1-11eb-1c50-634ba7a9a7bb
# ╠═f3fce35a-7b77-11eb-0935-13399babee72
# ╠═f80d6c4e-7b77-11eb-31bc-d55c9905b0a6
# ╟─fbf35d0a-7b77-11eb-3b36-194f011f00ab
# ╟─0b370ade-7b5d-11eb-14ae-8fb84aa69edc
# ╠═19fd463c-7b5d-11eb-342c-293f27a0f396
# ╠═2f54327a-7b5d-11eb-0352-e713aba2460d
# ╟─3a62156a-7b5d-11eb-1b66-a9d7d71d1051
# ╟─96d3fed8-7b7b-11eb-1434-07da91c27aca
# ╠═c5f3dde8-7b7b-11eb-3a14-bbe7f374c856
# ╟─2a044f08-7b7b-11eb-19fc-9fe29ad589bf
# ╠═8b29e5b4-7b7f-11eb-2679-bd248ed2d88b
# ╟─0f9373f0-78c5-11eb-26ba-05bea8874141
# ╠═3933c4ee-78c5-11eb-2692-7989afac7ebb
# ╟─6dda76f6-78c6-11eb-0fe0-a36a5663dfd7
# ╠═add65ba2-78c7-11eb-0eb6-771466c0c417
# ╠═d482635e-7ab4-11eb-0d53-398ca84dab54
# ╠═70ab6652-7ab4-11eb-1309-8f833821ef88
# ╟─3cf74998-78c5-11eb-1191-a9f56dbadda6
# ╠═3cfda932-78c5-11eb-3ad9-512459ffebc4
# ╟─3d11b968-78c5-11eb-3962-f50b8e05803d
# ╟─3d198bac-78c5-11eb-35af-87b99956e4de
# ╟─9d74de0e-7ac6-11eb-10c2-e3ac9cf17111
# ╟─7bd38294-7b7c-11eb-0329-99e46b366e1b
# ╟─64484bfa-7b7c-11eb-1d40-296fcf167f53
# ╟─0ff8fb0c-7b84-11eb-2a6a-0dedd44b3c7b
# ╟─cc6069f2-7b7c-11eb-3cb9-afbd773708e2
# ╟─18b49d34-7b80-11eb-3b12-3774bdf47ae6
# ╟─b2b0b842-7b57-11eb-1727-95c07b5b2de6
# ╟─7779ed4c-7b7d-11eb-034b-5fe2925df5df
