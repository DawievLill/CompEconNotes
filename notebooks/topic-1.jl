### A Pluto.jl notebook ###
# v0.12.21

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
md"# Floating point arithmetic

The goal of this session is to provide a brief overview of some of the elementary ideas used in numerical methods and approximation procedures. The discussion in this section mostly focusses on sources of inexactness in computational science, which are referred to as approximations. We are specifically referring to approximations of real numbers with a finite number of digits by means of floating point numbers. 

Errors generally arise from two basic aspects of numerical operations, rounding and truncation. We will talk about truncation errors at another time, for now the focus is on rounding errors. In order to appreciate rounding errors, we need to better understand floating point numbers and computer arithmetic. "

# ╔═╡ e6a2f610-7ac6-11eb-2920-0ba6b5af0ee6
md" ## Floating point numbers "

# ╔═╡ f1f991e0-7864-11eb-1b11-abf2aa055858
md" If you did not have time to go through any of the material on Julia so far, I recommend that you look [here](https://computationalthinking.mit.edu/Spring21/basic_syntax/) just to get a basic idea of the syntax. 

In this section we introduce the idea of floating point numbers. It is not crucial for the rest of the work, but it is useful to understand. Entire books have been written on floating point numbers, so this is just a brief introduction. Here is a nice video that motivates the usage of floating point numbers. I like the informal way in which it is discussed in the video."

# ╔═╡ e25b1a6a-7b75-11eb-06de-bdfec783caa1
html"""

<iframe width="680" height="400" src="https://www.youtube.com/embed/PZRI1IfStY0" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

"""


# ╔═╡ e6015b5c-7b75-11eb-2671-b19f7e366107
md"You might have heard the term `binary` before. This refers to a base 2 number system. Each digit can take on the value of 0 or 1. 

`bit` stands for `binary` + `digit`. 

A `byte` is equal to 8 `bits`, a `kilobyte` is $$10^3$$ `bytes` and a `megabyte` is $$10^6$$ `bytes`. This pattern continues for gigabytes, terabytes, petabytes, exabytes, zetabytes and so forth."

# ╔═╡ 85552d00-7865-11eb-1a94-3d5d52ed1772
md" As an aside, we can determine the amount of memory used by an object with the `Base.summarysize` command in _Julia_. I am not aware a specific function that does this in _Matlab_, but there is a `whos()` funtion that returns the size of all variables in the workspace. The _Julia_ equivalent of `whos()` is `varinfo()`."

# ╔═╡ 84df01a6-7861-11eb-2e3d-bd2e48f379f0
xx = rand(100, 100);

# ╔═╡ 6db2ed2a-7862-11eb-0053-99ca2c0cc1e6
Base.summarysize(xx);

# ╔═╡ 387561c8-7ce9-11eb-3458-03ad24bc43c3
md" ### Why do we need floating point numbers?"

# ╔═╡ 305a2bfe-7863-11eb-06c7-39c0544baece
md"A fixed point number system is a computer model for integers $$\mathbb{Z}$$, while the floating number system is a computer model for real numbers $$\mathbb{R}$$. 

We want to use a finite set of floating point numbers, denoted by $\mathbb{F}$, to represent $\mathbb{R}$. 

The fact that the real numbers are infinite along two dimensions (continuous and uncountable) is problematic in computing. Our objective is to understand the set $\mathbb{F}$, how rounding occurs when you operate on floating point values and also how rounding errors accumulate, and how you analyse the accuracy of numerical algorithms."

# ╔═╡ 6444bf22-7cbb-11eb-045b-7d7c90a508cb
md" Let us take a quick look at how the set $\mathbb{F}$ is represented on a computer with the following slider. We start with integers and then move on to floating point numbers.  "

# ╔═╡ b958066a-7c69-11eb-32d8-cd9430680e78
@bind m Slider(1:100)

# ╔═╡ c7b5a064-7c69-11eb-229e-3fc32f8632d4
m

# ╔═╡ 68335aa4-7c6a-11eb-3492-c3e8523e37a8
md" The `bistring()` command gives us an idea of the number is represented in the computer. Integers are represented in basic binary. Check what happens to the last few digits when you move the slider. "

# ╔═╡ 5e84490a-7c6a-11eb-1047-17f2a87d2ad8
bitstring(m)

# ╔═╡ a78a7548-7c6a-11eb-2246-1dafe0c2d816
md" Now we convert this integer to a floating point number in order to see its representation in the computer. In our minds these are basically the same, but we will soon see that these different types are seen differently by computers. You will also notice that the bitstring representation for floats and integers are not the same if you play around with the slider. "

# ╔═╡ d5a8db32-7c69-11eb-009e-0992320b2d2a
float(m)

# ╔═╡ e31c82d2-7c69-11eb-2608-9fadcd85a669
bitstring(float(m))

# ╔═╡ ecec8a9c-7c72-11eb-307c-a3eccc978923
md" Here are some other examples of floating point numbers. "

# ╔═╡ 045574b2-7c73-11eb-255f-0b43e4024c72
1.5e7  # Floating point value of 1.5 * 10^7

# ╔═╡ 0ada5afe-7c73-11eb-3f48-6ff9cc36362b
r = 1 / 49  # Division of two integers produce a floating point value

# ╔═╡ 36534b10-7c73-11eb-07a5-37a4fb5cc343
md" Note that $1 / 49 \notin \mathbb{F}$, our value given by $r$ is a rounded version of $1/49$. If we multiply by $49$ we see that it gives something close to but not exactly equal to $1$."

# ╔═╡ e720262c-7c73-11eb-0838-ab15f18a90dc
r * 49

# ╔═╡ edf67168-7c73-11eb-2cd5-d73de1a10454
1 - r * 49

# ╔═╡ 0a9dea26-7c74-11eb-037b-b7d58b25cd03
md" This difference is about $10^{-16} \approx 2^{-52}$ since the default precision in _Julia_ is **double precision**. " 

# ╔═╡ 5bdb0f86-84ee-11eb-2869-d9f08e860d61
bitstring(r)

# ╔═╡ cafaabec-84ee-11eb-35f5-15a895dbe372
md" Ok, so what causes this? Well, how are numbers represented in a computer: A bit string of 64 characters. 

Not simple binary, but a code where the position a digit occupies means something very precise. The first digit contains the `sign` of the number (positive or negative), the next 11 digits contain the order of magnitude (the `exponent`), and the final 52 digits, the fractional value in mathematical notation (in base 2/binary) that represents the number (the `mantissa`). 

This allows an extremely flexible way of (approximately) representing an enormous (but finite) number of real numbers exactly. The problems arise with those real numbers that cannot be represented exactly."

# ╔═╡ c2ebed22-84ed-11eb-0a10-a135b1823f78
md"The following image on double floating point format helps visualise the idea of a 64-bit precision representation as described. The first bit is a `sign` bit. Then there is the 11 bit `exponent` field and finally the 52 bit `mantissa`. Compare against the bit string for `r`."

# ╔═╡ f6c2e2b8-84ed-11eb-3026-890a1d21c649
md" ![floating point numbers](https://upload.wikimedia.org/wikipedia/commons/thumb/a/a9/IEEE_754_Double_Floating_Point_Format.svg/1024px-IEEE_754_Double_Floating_Point_Format.svg.png)"

# ╔═╡ 45d983b6-84ee-11eb-327c-55f5a3ca3341
md"We will continue this discussion on double precision and floating-point number representation after a quick detour into a discussion on types."

# ╔═╡ d66c205c-7c54-11eb-1500-2deb82c42c6c
md" #### Quick detour into types "

# ╔═╡ d5191c50-7b55-11eb-20c9-896d13cd253e
md"Before we move on with floating point numbers, let us quickly talk about types. In _Julia_ (and most other programming languages) `1` and `1.0` are different values, because they have different types. One could also think of other ways to represent one in different ways (all of which are equally valid). Look at the following possible representations of one."

# ╔═╡ de731d1a-7c53-11eb-10d6-3572e8e7ee70
one = [
	1,
	1.0,
	"one",
	1//1,
	[1 0; 0 1],
]

# ╔═╡ e3baefd2-7c53-11eb-25b2-23f4bb7e6908
md" **Types** are incredibly important and we will talk about them frequently during this reading group. Each of the items in this list is a specialised representation of one. The diffence between these representations are clear to me, but how does a computer see this? We can determine the type of any object by using the `typeof()` command."

# ╔═╡ ee60791a-7b55-11eb-0bed-2d54e6d2f18b
typeof(1)

# ╔═╡ f92a2226-7b55-11eb-2dd1-bbdb6d228fef
typeof(1.0)

# ╔═╡ 74d937c6-7c54-11eb-093e-a3c49b0340c1
md" Instead of doing this `typeof()` command for each value in the list, we can `broadcast` the `typeof()` operation across all values using the `.` operator as seen below."

# ╔═╡ 689c2086-7c54-11eb-2e4d-f7b6c988d7c3
computer_ones = typeof.(one)

# ╔═╡ 4140fd08-7c55-11eb-1583-7d3c3de7c8ef
md" We can clearly see that the computer believes that these are all different types "

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
md"Double precision numbers are the dominant data type in scientific computing. `Float64` is the type for double precision numbers in _Julia_. As far as I know, _Matlab_ only uses the double precision floating point number. These numbers utilise $64$ bits of memory."

# ╔═╡ 237a5fe0-7c6b-11eb-2fa1-7f3135c6faa0
md" ### Floating point representation " 

# ╔═╡ 56547ac2-7cba-11eb-2698-a13ef1195c67
md" This first part is perhaps a bit more techical and you can just skim through on your first reading. The objective is to help you understand how numbers are represented. In a computer the real numbers are represented by a floating-point number system with a fixed number of digits. 

The floating-point number system is similar to scientific notation and $\mathbb{F}$ is characterised by four integers, namely the base $b$, precision $p$ and exponent range $[L, U]$. 

Any floating point number $x \in \mathbb{F}$ has the form 

$$x = \pm \left(\frac{d_1}{b} + \frac{d_2}{b^2} + \cdots + \frac{d_{p}}{b^{p}}\right) \times b^{e} = \pm (r) \times b^{e}$$ 

Here $d_i$ is an integer such that $0 \leq d_i \leq b$ with $i = 1, \ldots, p$.

In this representation $e$ is an integer such that $ L \leq e \leq U$. 

The part in parentheses, which is represented by a string of $p$ base-$b$ digits ($. d_1\ldots d_{p}$), is called the `significand` or `mantissa` (which we often indicate with an $r$). You should note that $1 \leq r < b$.

The `sign` is represented by $s$. The number of digits of $p$ is the `precision` and $e$ is the `exponent`. 

As an example, if the base is $10$, then we have that:

$x = \pm (r) \times 10^{e}$

with $1 \leq r < 10$. If we want to represent $2347$, then we have $2.347 \times 10^{3}$."

# ╔═╡ 735dbd02-7f85-11eb-197a-9d5bd5750b35
md" ##### Double precision"

# ╔═╡ 2f56aafa-7f86-11eb-0e4c-51b5754a86ec
md"Floating-point numbers are traditionally represented with base of $2$, which means that we have:

$x = \pm (r) \times 2^{e} = s (r) \times 2^{e}$

with $1 \leq r < 2$. We replace the $\pm$ with $s$ to indicate the `sign`. "

# ╔═╡ 8f4efeec-7f85-11eb-20a5-d705eb5d9b2e
md" Double-precision binary floating-point is a commonly used format on PCs, due to its wider range over single-precision floating point, in spite of its performance and bandwidth cost.  "

# ╔═╡ f761979e-835a-11eb-364b-09fd2e90b6d8
md"For the binary version of the number system we that $r$ is a **normalised** `mantissa`. This means that $r$ will have the following form:

$r = 1.\text{(fractional part)}$

Normalisation forces the integer part of the mantissa to be exactly one and therefore is no reason to store it in the computer. The fractional part can be whatever we like it to be. "

# ╔═╡ c2f91cf4-835d-11eb-3b5a-fb7f329f4697
function dec2bin(x::String)
    bx = parse(BigFloat, x)
    xneg = bx >= 0 ? false : true
    bx = BigFloat(xneg ? -bx : bx)
    pre, post, n = "", "", div(nextpow(2, bx), 2)
    while bx > 0
        a, bx = divrem(bx, BigFloat(n))
        if n >= 1
            pre *= a > 0 ? '1' : '0'
        else
            post *= a > 0 ? '1' : '0'
        end
        n /= 2.0
    end
    (xneg ? "-" : "") * pre * "." * post
end

# ╔═╡ d4482e58-835d-11eb-1cf1-b1a6d9c5c0b2
dec2bin("13.25")

# ╔═╡ e1e354b0-835b-11eb-0682-9da5144f8f99
md" Let us illustrate this quickly with an example. If you take the number $13.25$ in binary you will see that it gives $1101.01$. Look at the code above for an illustration. " 

# ╔═╡ 0859b1c8-835e-11eb-3be5-2163e163e365
md"The integer part of this binary number is $1101$ and $01$ is the fractional part. One could represent $13.25$ as $1101.01 * (2^0)$, but this is not normalised since the integer part is not $1$. However, we can shift the `mantissa` to the right one digit if we increase the exponent by one:

$ 1101.01 * (2^0) = 110.101 * (2^1) = 11.0101 * (2^2) = 1.10101 * (2^3)$.

This representation in the last equality is the **normalised** form of $13.25$. "

# ╔═╡ 46c6770e-831b-11eb-067e-21334c077e86
md" Continuing from our depiction in scientific notation in the previous section, we can now use the binary base to obtain the representation for the double precision number as the following:

$s (r) \times 2^{e} = s \left(1.d_2\ldots d_{p+1}\right)_{2} \times 2^{e} = s \left(1 + m/2^{52}\right)\times 2^{e}$

where $m = {d_1}\times{2}^{-1} + {d_2}\times{2}^{-2} + {d_3}\times{2}^{-3} \cdots$"

# ╔═╡ f8ca7706-831c-11eb-1ac3-77fa1483b53f
md" Let us try and transform this representation so that it reflects the IEEE 64-bit floating-point representation. We start with the sign $s$, which can be $\pm 1$. Another way to write our representation:

$(-1)^s \left(1.d_2\ldots d_{53}\right)\times 2^{e}$" 

# ╔═╡ c5b1fb6a-8327-11eb-2e9a-4fc8d2b3eac2
md"We know that $s$ needs to be represented in binary, since our base is $2$. In this representation we have that if $s = 0$ then the number will be positive and if $s = 1$, then the number is negative."

# ╔═╡ 58110390-8329-11eb-1c11-8b24865d91ff
md" Finally we have to deal with the exponent. The main problem with the current representation is that negative exponents could pose a problem in comparisons. Let us look at an example to see why this is the case."

# ╔═╡ 1c178a14-831d-11eb-2380-93561ebc9bd8
bitstring(1.0 * 2^(-1))

# ╔═╡ dbb19784-831d-11eb-105e-ed8117a41c6b
bitstring(1.0 * 2^(1))

# ╔═╡ ed24d814-831d-11eb-1018-1d0c85976c0b
md" With this representation of the numbers $1.0 * 2^{-1}$ and $1.0 * 2^{1}$ we see that the first exponent shows a larger binary number. The exponent here is represented by digits 2 though 12 in the bitstring. This makes direct comparison difficult."

# ╔═╡ 260bba26-831e-11eb-14b1-47d4ae08397d
md" This motivates the usage of **biased notation** for exponents." 

# ╔═╡ 4107719c-831e-11eb-1766-cb889deb02de
md" If the real exponent of a number is $X$ then it is represented as $X$ + `bias`. The bias in IEEE double precision is $1023$, therefore the exponent $-1$ is represented as $-1 + 1023 = 1022$. This means that the actual exponent is found by substracting the bias from the stored exponent."

# ╔═╡ c6a576d2-7f85-11eb-1e11-43906c0c1187
md" Now we have our final IEEE represenation with a given **biased** `exponent` $e$ and a 52-bit `mantissa`:

$(-1)^\text{sign} \left(1.d_{2}d_{3}\ldots d_{53}\right)\times 2^{e-1023}$

or alternatively, 

$(-1)^\text{sign} \left(1 + m/2^{52}\right)\times 2^{e-1023}$"

# ╔═╡ 531698fe-7f85-11eb-0f2d-bbc361c47c72
md"Remember the image from before with respect to floating point numbers. The first bit is a `sign` bit. Then there is the $11$ bit `exponent` field, which is an unsigned integer from $e = (00\ldots 0)_2 = 0$ to $e = (11\ldots 1)_2 = 2047$."

# ╔═╡ ee13edce-84ed-11eb-394d-85367178ceec
md" ![floating point numbers](https://upload.wikimedia.org/wikipedia/commons/thumb/a/a9/IEEE_754_Double_Floating_Point_Format.svg/1024px-IEEE_754_Double_Floating_Point_Format.svg.png)"

# ╔═╡ 4fb1b54a-8320-11eb-0cd4-b38d5e72b1b1
md"We use the exponent in biased form as explained before. The biased form means that an exponent value of $1023$ represents the actual value of zero. One calculates the value of the exponent by subtracting the bias for `exponent` to get an exponent value in the range $[−1023, 1024]$. The bounds for the exponent are: $0 \leq e \leq 2047$. The final part is the $p = 52$ bit `significand`."

# ╔═╡ a54ba598-831a-11eb-2dcf-dbc008024c55
md" With double precision we can then see that the floating-point number occupies 64-bits. There is a compromise between the size of the `significand` and the `exponent`. The chosen size for the exponent provides a largest positive number of $2^{1024} \approx 2 \times 10^{308}$ and a smallest positive number of $2^{-1022} \approx 10^{-308}$. In the case that the exponent is too large to be represented in the exponent field we have **overflow**. If the number is too small to be represented in this exponent field we have **underflow**. "

# ╔═╡ df023fb6-7f70-11eb-215a-7d834c2e7b56
md" ##### Looking under the hood "

# ╔═╡ 06bfbdd0-7c6a-11eb-033b-4df5aaeb3220
md"""
The bit of code that we have below will give us the different components of this IEEE 64 bit representation (double precision).
"""

# ╔═╡ fbe88cde-7c69-11eb-2c11-f744889098a8
ieee0(x::Float64) = [bitstring(x)[1:1], bitstring(x)[2:12], bitstring(x)[13:64]]

# ╔═╡ 9ba79240-7c6c-11eb-3385-7917a37df361
md" Remember that with floating point numbers we are working with a computer representation of real numbers. One of the properties of $\mathbb{R}$ is that it is continuous. Let's focus our attention on a range of values in the interval $[1, 2)$. This will be the interval if our value for the `exponent` is set to zero. "

# ╔═╡ 74bec6fa-7c6d-11eb-27d0-c972ad74515c
@bind n Slider(1: 2^(52), show_value = true)

# ╔═╡ 803ff050-7c6d-11eb-1bc3-ff5335e4b361
md" Generally the spacing between elements of $\mathbb{F}$ in $[2^{e}, 2^{e+1})$ is $2^{e-p}$. Therefore you will see that the values for the slider range from $1:2^{52}$. This slider provides us with a range of values between 1 and 2. For the real numbers there are an infinite amount of values between 1 and 2, but we can only present a limited amount on a computer. In fact, $2^{52}$ is maximum that we can present." 

# ╔═╡ 14591b22-7c6e-11eb-008f-27a2ca96c8e9
q = 1 + n/2^52

# ╔═╡ 29b3d9bc-7c6e-11eb-2702-53b1dddb63ea
md" In this case $q$ represents all the values between 1 and 2 that the computer can take. In the code below we see the different components of the floating point number for values in our selected interval."

# ╔═╡ 2111ff28-7c6e-11eb-1ff1-f5017b3889fa
ieee0(q)

# ╔═╡ 955462da-7c6f-11eb-0fe9-a7f3913bf97d
md" The problem for us is that these are not completely readable, since we are used to working with base-10. So let's parse this information using another bit of code." 

# ╔═╡ bfaf9f4c-7c6f-11eb-01b3-8334e46772b6
ieee(x) = [ parse.(Int, ieee0(x), base=2) ieee0(x)]

# ╔═╡ 988d1952-80fb-11eb-3f37-25d71c811413
md" So before we look at the output of the function let us remember how to read this:

$(-1)^s  * 2^{(e-1023)} * (1 + m/2^{52})$
"

# ╔═╡ d0a53b18-7c6f-11eb-0afd-0faa5c85a613
ieee(q)

# ╔═╡ 5f54346a-7866-11eb-3df4-658de2c2ce5f
md"### Floating point arithmetic"

# ╔═╡ 600b116c-84eb-11eb-00e0-f14d7b241d24
md" Here is a nice video that explains most of the concepts on floating-point representation that we covered and some discussion on computer arithmetic (although it talks about algebra). It is a bit more technical, but it provides some interesting insights."

# ╔═╡ 1d08b996-84eb-11eb-0736-510e5f7cdd1f
html"""

<iframe width="680" height="400" src="https://www.youtube.com/embed/pQs_wx8eoQ8" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

"""

# ╔═╡ 5cf34e4e-7867-11eb-18ff-01b71b0d5934
md"The basic issue is that, for computer arithmetic to be fast, it has to be done in hardware, operating on numbers stored in a fixed, finite number of digits (bits). As a consequence, only a finite subset of the real numbers can be represented, and the question becomes which subset to store, how arithmetic on this subset is defined, and how to analyze the errors compared to theoretical exact arithmetic on real numbers.

Let us work with some examples to show what problems we might encounter. "

# ╔═╡ 7c389718-8321-11eb-2e05-d3969090fe0e
md" #### Rounding error"

# ╔═╡ 0926daac-7ccb-11eb-21ca-391422834af3
md" The results of the following floating point calculations are often surprising to people. "

# ╔═╡ bd2cb4c8-8323-11eb-34dc-a1cdf2cccb25
(sqrt(3))^2

# ╔═╡ d12e4e0a-8323-11eb-0d74-6b574609a909
md"Shouldn't the answer be 3? Here we have the problem that in the original process of taking the square root, we are rounding off the answer. Once we try to square this number we are taking the square of a number that has been rounded. " 

# ╔═╡ c54e39f0-7ccb-11eb-3d69-331302b5def2
5/6

# ╔═╡ c950647e-7ccb-11eb-028f-d1c7c065342f
md" Shouldn't the last digit be a 3? The reason for the trailing 4 is because this is the 64-bit float that is closest to 5/6 when converted back to decimal. By the way, in Matlab you won't notice this 4 since it only prints the first 15 digits. Matlab is sneaky that way!"

# ╔═╡ c89d83a4-7ccb-11eb-210e-355010a04b45
2.6 - 0.7 - 1.9

# ╔═╡ e17b1136-7ccb-11eb-1c45-e9de3d9cd831
md" Shouldn't this answer be 0? In this case none of the three decimal numbers has an exact 64-bit binary representation, so each must be rounded to the nearest 64-bit float before doing arithmetic. The number that we see in the output is actually $10^{-16} \approx 2^{-52}$."

# ╔═╡ 11ea91ce-7cd0-11eb-2aa5-3df9a4b24db5
md" These are not bugs in _Julia_ and you will encounter them in all programming languages that adopt the IEEE-standard 64-bit binary representation of floating point numbers. " 

# ╔═╡ 479d8d1c-8324-11eb-2d03-8dfbf74ea283
md" Conducting computer arithmetic differs fundamentally from the way that you would do it by hand. There are two obvious reasons for this. 

1. Only finite many numbers can be represented in a computer.
2. A number represented in a computer can only have finitely many digits."

# ╔═╡ 8fdc5676-8324-11eb-01c6-b1f82b890a40
md" Let us look at some examples of computer arithmetic to make the problems clear."

# ╔═╡ ac567efc-7867-11eb-07c3-21989feec15c
begin
	a = 0.1
	b = 0.1
	c = 0.1
	
	a + b + c == 0.3
end

# ╔═╡ 31795f9a-7b7a-11eb-1786-89aad77aff4b
md" As we mentioned in the example before, this problem is due to rounding. Rounding is necesarry when a number has more than $p$ `significand` bits. _Julia_ offers several rounding modes, with the default being `RoundNearest`. In our example here, the number `0.1` in the decimal system cannot be represented accurately as a floating point number. " 

# ╔═╡ de5fc23a-7b7a-11eb-12b0-a513044b39a6
ieee(0.1)

# ╔═╡ e7475a98-7b7a-11eb-00a2-63e52e403c16
md"Do you see the repeating pattern in the `significand`? Even thought `0.1` can easily be represented with only two base-10 digits, it requires an infinite number of binary digits, which is cut off. This quotient of two integers has a nonterminating binary expansion is not a binary floating-point number. This is a confusing aspect of floating-point arithmetic. We have a real arithmetic operation on two floating-point numbers that does not result in aother floating-point number. 

If a number that is not representable as a floating point number is entered into the computer then it must be rounded to obtain a floating-point number. It is possible to then have rounding errors on human-readable decimal values for both input and output. There is something called [decimal floating point](https://en.wikipedia.org/wiki/Decimal_floating_point) that avoids this particular issue, but it is slow and only used for relatively specialised purposes. We will cover this in another example toward the end of the session. One could also use **BigFloats** that have an arbitrary amount of precision. "

# ╔═╡ e54c99d8-7c74-11eb-29cc-f5d7f64c4937
md" Given the previous example why do you think the following example is different?"

# ╔═╡ 0977c012-7c75-11eb-3cf4-e9a1bfca6541
1.0 + 1.0 == 2.0

# ╔═╡ 314c5f4e-7c75-11eb-2b48-5787c4c781c9
md" This is since basic mathematical operations like addition and subtraction are guaranteed to be true for integer arithmetic in floating point until you exceed the largest representable integer in your precision. This is a more subtle point and requires some further reading in the fundamental axioms of floating point arithmetic."

# ╔═╡ 3310e738-7ac1-11eb-38ea-2f1cfc2f0090
md"Let us illustrate with another example that exact arithmetic and computer arithmetic don't always give the same answers. Consider the following:"

# ╔═╡ 8f365958-7ac1-11eb-145b-69bff762b7a8
begin 
	e = (1e-20 + 1) - 1
	f = 1e-20 + (1 - 1)
	e == f
end

# ╔═╡ 9fb90230-7ac1-11eb-1c50-634ba7a9a7bb
md" The value we obtain for `e` is actually incorrect."

# ╔═╡ f3fce35a-7b77-11eb-0935-13399babee72
e

# ╔═╡ f80d6c4e-7b77-11eb-31bc-d55c9905b0a6
f

# ╔═╡ fbf35d0a-7b77-11eb-3b36-194f011f00ab
md" Even though the statements are mathematically the same (by the associative property of addition and subtraction), if we compare the values we see that `e` > `f`. In this case this is because adding numbers of different magnitudes does not always work like you would want. In the example, when we added $10^{-20}$ to $1$, it got rounded away. This means the floating point arithmetic is not associative in general. "

# ╔═╡ ad6d4176-8321-11eb-2cf0-89aa37f483f5
md" #### Catastrophic cancellation"

# ╔═╡ 0b370ade-7b5d-11eb-14ae-8fb84aa69edc
md" Unfortunately rounding errors are not the only cause for concern in floating-point arithmetic. Another class of problems using floating point arithmentic is illustrated in the following example:"

# ╔═╡ 19fd463c-7b5d-11eb-342c-293f27a0f396
begin 
	hh = 100000.2
	ii = 100000.1
end;

# ╔═╡ 2f54327a-7b5d-11eb-0352-e713aba2460d
hh - ii

# ╔═╡ 3a62156a-7b5d-11eb-1b66-a9d7d71d1051
md"Think about what you would expect the answer to the calculation to be. Did you get it right?"

# ╔═╡ 96d3fed8-7b7b-11eb-1434-07da91c27aca
md" What about the following example. Think about what the result should be and what the computer calculates. Did you get this one right? "

# ╔═╡ c5f3dde8-7b7b-11eb-3a14-bbe7f374c856
begin
	j = 2.0^30 # Very large number
	k = 2.0^-30 # Very small number
	j - k == j
end

# ╔═╡ 2a044f08-7b7b-11eb-19fc-9fe29ad589bf
md"Both of the examples above showcase the idea of **catastrophic cancellation**. In the first instance the subtraction of two nearly equal numbers eleminates significant digits (since the numbers can't be precisely represented by the machine) and in the second case there is a potential loss of precision when subtracting distant real numbers that are represented by floats.

This type of error is a common source of huge floating point errors."

# ╔═╡ 8b29e5b4-7b7f-11eb-2679-bd248ed2d88b
isapprox(100000.2 - 100000.1, 0.1)

# ╔═╡ bc86d990-83d2-11eb-1401-e1ac3ac6b36c
md"For the first example we see that the values are approximately equal, as shown above."

# ╔═╡ 95750e6c-83d2-11eb-31e4-673d932bb999
md" ##### Calculating returns on investment "

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
end

# ╔═╡ 6dda76f6-78c6-11eb-0fe0-a36a5663dfd7
md"We have done the calculation using floating point numbers, now we use the precise decimal representation for a comparison. Using this decimal representation is slow, so we would rather like to use floating point arithmetic. However, for the purpose of the example let us continue. "

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
end

# ╔═╡ 3cf74998-78c5-11eb-1191-a9f56dbadda6
md"Compare the results and see if there is a significant difference. "

# ╔═╡ 3cfda932-78c5-11eb-3ad9-512459ffebc4
diff = sum_1 - Dec64(sum)

# ╔═╡ 3d11b968-78c5-11eb-3962-f50b8e05803d
md" Do you expect there to be any difference between the answers?"

# ╔═╡ 6ff7cb16-83d2-11eb-0155-e38cae990886
md"##### Comparing functions"

# ╔═╡ 441ea182-8361-11eb-2fef-11c314784df1
md" In this next example we will make use of two functions to see where the loss of significant digits can be really worrisome. Consider computing the function:

$p(x) = (x-1)^5$

One can expand the polynomial to find an *equivalent* function:

$q(x) = x^5 - 5x^4 + 10x^3 - 10x^2 + 5x - 1$

However, we run into problems when $x \approx 1$. The alternating signs lead to subtractions of two very close numbers, which leads to large error!"

# ╔═╡ fd34d75e-8363-11eb-1ddd-efaf1569a22b
md" Our code below generates plots for $p(x)$ and $q(x)$ on the interval $1-10^{-10} \leq x \leq 1 + 10^{-10}$."

# ╔═╡ e688f190-8361-11eb-1fe8-2bc7c5295926
X = 1.0 - 10^(-10) : 10^(-12) : 1.0 + 10^(-10)

# ╔═╡ 4a4ab3ca-8363-11eb-1c65-cb866fc822ee
md"**Note:** You can use the `@` symbol to invoke `macros`. In this case the macro allows me to broadcast the `.` operator across all arithmetic operations. This means I don't repeatedly have to type `.` before each arithmetic operation."

# ╔═╡ 53b3db4a-8362-11eb-1948-7b49d75d9f99
P = @. (X - 1)^5

# ╔═╡ b77ec01c-8363-11eb-037f-833825146fa5
plot(P,
    xlabel = "p(x)",
    xlims = (0.9999,1.001)
)

# ╔═╡ 75d0456a-8362-11eb-244b-7f53c80e5bc8
Q = @. X^5 - (5 * X^4) + (10 * X^3) - (10 * X^2) + (5 * X) - 1 

# ╔═╡ a7de2550-8363-11eb-038d-696d2e53f139
plot(Q,
    xlabel = "q(x)",
    xlims = (0.99999,50)
)

# ╔═╡ caeae7de-8363-11eb-003d-83d7fbea4217
md"One would expect that the result would be pretty flat in the region around the value of $1$, but in the graph of $q(x)$ we see this strange behaviour. There is large error here with wild oscillations, which does not capture the properties of the function at all. "

# ╔═╡ c8daad14-8361-11eb-169d-e36ccbd815ba
md" ##### Not so catastrophic cancellation?"

# ╔═╡ 9fb8ee4a-7cc5-11eb-0a93-2fea93a38c75
md" Catastrophic cancellation is not inevitable! We can often avoid it by simply re-arranging the calculation. Let us look at a final example to see how to do this. "

# ╔═╡ 5974daa4-7cc6-11eb-097a-1303b35124dc
md" Suppose you want to calculate the function $e^x - 1$ using floating point arithmetic. When $|x| \ll 1$, we have that  $e^x \approx 1$ and so the calculation $e^{x} - 1$ will experienc catastrophic calculation."

# ╔═╡ c8e07448-7cc6-11eb-1e2d-15d48c74d921
x = 2.0^-60

# ╔═╡ d2df623a-7cc6-11eb-2dfd-c57bb3dfddcd
exp(x)

# ╔═╡ df91c5f2-7cc6-11eb-17c4-f9129b2fc701
exp(x) - 1 # naive algorithm leads to catastrophic cancellation

# ╔═╡ ea286192-7cc6-11eb-3274-6fa2ee43cc44
md" The correct answer is in fact the following "  

# ╔═╡ f94f0b80-7cc6-11eb-1ad6-c134a7ae5590
Float64(exp(big(x)) - 1)

# ╔═╡ 03468256-7cc7-11eb-3b5f-a16deac85f21
md" You can also see this using the Taylor expansion of $e^{x}$, which we can use to calculate this function accurately for small $x$"

# ╔═╡ 439e7f0e-7cc7-11eb-0080-c970dd4196ba
x + x^2/2 + x^3/6 # 3 terms is more than enough for x ≈ 8.7e-19

# ╔═╡ 4efe834c-7cc7-11eb-1060-9137f0a8da19
md" The key is to rearrange the calculation to perform the cancellation analytically, and only use floating-point arithmetic after this is accomplished. In fact, Julia's standard library (and scientific-computing libraries in other languages) provides a function called `expm1(x)` that computes $e^{x}−1$ accurately for all values of x. "

# ╔═╡ 64337256-7cc7-11eb-2ef3-4304c425d51d
expm1(x)

# ╔═╡ 3d198bac-78c5-11eb-35af-87b99956e4de
md"We have to realise that real numbers are only presented with a certain precision. Most of the time we don't need to worry, but we want to make sure our code is robust. In other words, we shouldn't be surprised by the outcome of our code. The **numerical stability** of our code is important! "

# ╔═╡ 000779f6-7f6c-11eb-012b-6b67b412eb5b
md" ## Error analysis"

# ╔═╡ 1029caaa-7f6c-11eb-0121-5b4b86a9a832
md" One of the topics we did not have time for is the study of the effects of approximations on the accuracy and stability of numerical algorithms. This is traditionally referred to as error analysis. We will not undertake a deep look into error analysis in this section, but it it something that might appear in other sessions, so let us introduce some basic concepts."

# ╔═╡ 5f093552-7f6c-11eb-3ea4-5b0032275023
md" Two types of error are commonly used. Absolute and relative error, which are defined as:"

# ╔═╡ 7a843716-7f6c-11eb-3656-a1014b9984ad
md" `Absolute error` = `approx value` - `true value`. "

# ╔═╡ 92dda336-7f6c-11eb-2e80-d959593aff78
md" `Relative error` = `absolute error` / `true value`. "

# ╔═╡ b8c2920a-7f6c-11eb-37fc-9fa12702d8f2
md" A completely erroneous approximation would correspond to a `relative error` of at least 1. This means that `absolute error` is greater than the `true value`. "

# ╔═╡ 4fb32b18-8330-11eb-1cc8-779e4a09332d
md" Consider a basic example. Recall the high-school approximation to the number π."

# ╔═╡ 63b924f0-8330-11eb-121a-c56a30993c81
p = 22/7

# ╔═╡ 90a198a8-8330-11eb-3b78-fff35e711256
md"Not all the digits displayed for $p$ are the same as those of π. As an approximation, its absolute and relative accuracy are:"

# ╔═╡ 6a25b7e0-8330-11eb-2258-c357b0033d47
abs_accuracy = abs(p-pi)

# ╔═╡ 7a9c88b0-8330-11eb-0938-39c0cca872ae
rel_accuracy = abs(p-pi)/pi

# ╔═╡ a92a2796-8330-11eb-210b-11fb46621064
md"We can also calculate the number of accurate digits in $p$"

# ╔═╡ 84e1f6d4-8330-11eb-3fc9-271654090046
accurate_digits = -log(10,rel_accuracy)

# ╔═╡ ef4606b8-8318-11eb-1bf5-5fee8d428bae
md" ## Machine epsilon"

# ╔═╡ 0698f0c8-8319-11eb-2035-191ff7afc846
md" This is another important concept, related to ideas of underflow and overflow, that I did not have time to explore. However, in future iterations we will perhaps talk about this."

# ╔═╡ 7bd38294-7b7c-11eb-0329-99e46b366e1b
md" ## Rules for numerical computation "

# ╔═╡ 64484bfa-7b7c-11eb-1d40-296fcf167f53
md" Some basic rules that we want to follow when doing computation.

1. Add small numbers together before adding larger ones
2. Add numbers of like magnitude together. This is called pairing.
3. Avoid subtraction of two numbers that are nearly equal."

# ╔═╡ 0ff8fb0c-7b84-11eb-2a6a-0dedd44b3c7b
md" Note that subtractive cancellation is one of the most common mechanisms introducing dramatic growth of errors in floating point computation "

# ╔═╡ cc6069f2-7b7c-11eb-3cb9-afbd773708e2
md" You could have saved a lot of time and just looked at these rules, but the examples give an idea as to why you might want to implement the rules in the first place."

# ╔═╡ 18b49d34-7b80-11eb-3b12-3774bdf47ae6
md" We will encounter more types of errors in approximation in the future, especially when we talk about calculating derivatives. Then we will introduce ways in which we can quantify the error with different metrics. However, I think this is enough information for one session. Your next step is to look at some of the exercises below."

# ╔═╡ b2b0b842-7b57-11eb-1727-95c07b5b2de6
md" ## Exercises 

For the first few exercises let's check whether floating-point numbers obey certain algebraic rules. For 2-5, one couterexample will suffice. All these initial problems are considered **[Easy]**. 

1. The associative rule for addition says `(x + y) + z == x + (y + z)`. Check this rule using `x = 0.1`, `y = 0.1` and `z = 1.0`. Explain what you find.
2. Do floating-point numbers obey the associative rule for multiplication: `(x * y) * z == x * (y * z)`?
3. Do floating-point numbers obey the distributive rule: `a * (x + y) == a * x + a * y`?
4. Is `0 * x == 0` true for all `x` (where `x` is a floating point number)?
5. Is `x / a == x * (1 / a)` always true?"

# ╔═╡ 9bd103ee-7f75-11eb-0f0d-b13097313f0a
md"6. **[Easy-Medium]**  In this question we will be looking at the standard formula for finding the roots of a quadratic equation, which is given by: 

$x = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}$

If some values of the coefficients in this formula we can run into problems. If the coefficients are very large or very small then the values of $b^2$ or $4ac$ may **overflow** or **underflow**. We did not deal with this issue specifically in this session, so you can go `Google` those ideas and read more on them. The basic idea is that we don't have good representations of values for infinity or zero. This question does not focus on those issues, rather we would like to think about another possible problem that might occur with the usage of this quadratic formula in floating-point arithmetic. 

6.1. How would one avoid cancellation between $-b$ and the square root?

6.2. Consider the following coefficient values to see why we might have a problem with our arithmetic. Select $a = 0.05010$, $ b = -98.78$ and $c = 5.015$. The correct roots, rounded to ten significant digits are: $1971.605916$ and $0.05077069387$. Compute the roots using the formula and see how they compare. What did you find?
"


# ╔═╡ b8783bca-7f75-11eb-070c-697f8d79c0e8
md"7. **[Medium]** Exercise still to be completed"

# ╔═╡ 7779ed4c-7b7d-11eb-034b-5fe2925df5df
md" 8. **[Hard]** In this problem we will use `for loops`, so be sure that you have an idea of how to use loops. We can expect floating point errors to accumulate randomly during computation, creating what is known as a **random walk**. 

On average we expect the errors to partially cancel out. Suppose you define a random sequence by $x_0 = 0$ and $x_n = x_{n-1} \pm 1$ with the signs chosen by tossing a fair coin for each $n$. et $\alpha_n$ and $\beta_n$ be the average value of $x_n$ and $|x_n|$ respectively over all such walks. Then a classic result of probability is that $\alpha_n = 0$ and $\lim_{n\rightarrow \infty} \frac{\pi \beta^{2}_{n}}{2n} = 1$. Perform a million random walks. Exercise still to be completed. "

# ╔═╡ Cell order:
# ╟─5b44b6ec-7863-11eb-1ed4-3d0c9eadd065
# ╟─66ca1514-7863-11eb-2fd9-c5a2e7d1a74d
# ╟─844fa292-7b48-11eb-0f9c-c521a4125ce6
# ╟─8598d0ca-7863-11eb-1549-21bc81a5cb1f
# ╟─e6a2f610-7ac6-11eb-2920-0ba6b5af0ee6
# ╟─f1f991e0-7864-11eb-1b11-abf2aa055858
# ╟─e25b1a6a-7b75-11eb-06de-bdfec783caa1
# ╟─e6015b5c-7b75-11eb-2671-b19f7e366107
# ╟─85552d00-7865-11eb-1a94-3d5d52ed1772
# ╠═84df01a6-7861-11eb-2e3d-bd2e48f379f0
# ╠═6db2ed2a-7862-11eb-0053-99ca2c0cc1e6
# ╟─387561c8-7ce9-11eb-3458-03ad24bc43c3
# ╟─305a2bfe-7863-11eb-06c7-39c0544baece
# ╟─6444bf22-7cbb-11eb-045b-7d7c90a508cb
# ╠═b958066a-7c69-11eb-32d8-cd9430680e78
# ╠═c7b5a064-7c69-11eb-229e-3fc32f8632d4
# ╟─68335aa4-7c6a-11eb-3492-c3e8523e37a8
# ╠═5e84490a-7c6a-11eb-1047-17f2a87d2ad8
# ╟─a78a7548-7c6a-11eb-2246-1dafe0c2d816
# ╠═d5a8db32-7c69-11eb-009e-0992320b2d2a
# ╠═e31c82d2-7c69-11eb-2608-9fadcd85a669
# ╟─ecec8a9c-7c72-11eb-307c-a3eccc978923
# ╠═045574b2-7c73-11eb-255f-0b43e4024c72
# ╠═0ada5afe-7c73-11eb-3f48-6ff9cc36362b
# ╟─36534b10-7c73-11eb-07a5-37a4fb5cc343
# ╠═e720262c-7c73-11eb-0838-ab15f18a90dc
# ╠═edf67168-7c73-11eb-2cd5-d73de1a10454
# ╟─0a9dea26-7c74-11eb-037b-b7d58b25cd03
# ╠═5bdb0f86-84ee-11eb-2869-d9f08e860d61
# ╟─cafaabec-84ee-11eb-35f5-15a895dbe372
# ╟─c2ebed22-84ed-11eb-0a10-a135b1823f78
# ╟─f6c2e2b8-84ed-11eb-3026-890a1d21c649
# ╟─45d983b6-84ee-11eb-327c-55f5a3ca3341
# ╟─d66c205c-7c54-11eb-1500-2deb82c42c6c
# ╟─d5191c50-7b55-11eb-20c9-896d13cd253e
# ╠═de731d1a-7c53-11eb-10d6-3572e8e7ee70
# ╟─e3baefd2-7c53-11eb-25b2-23f4bb7e6908
# ╠═ee60791a-7b55-11eb-0bed-2d54e6d2f18b
# ╠═f92a2226-7b55-11eb-2dd1-bbdb6d228fef
# ╟─74d937c6-7c54-11eb-093e-a3c49b0340c1
# ╠═689c2086-7c54-11eb-2e4d-f7b6c988d7c3
# ╟─4140fd08-7c55-11eb-1583-7d3c3de7c8ef
# ╟─a0592f30-7b78-11eb-2482-f5eaa3a66228
# ╠═6f5370c6-7b78-11eb-2076-45dc8e4908ae
# ╟─65b2d784-7b78-11eb-063b-1b7c30edd211
# ╟─4e406a14-7b56-11eb-109c-6136655fc096
# ╟─237a5fe0-7c6b-11eb-2fa1-7f3135c6faa0
# ╟─56547ac2-7cba-11eb-2698-a13ef1195c67
# ╟─735dbd02-7f85-11eb-197a-9d5bd5750b35
# ╟─2f56aafa-7f86-11eb-0e4c-51b5754a86ec
# ╟─8f4efeec-7f85-11eb-20a5-d705eb5d9b2e
# ╟─f761979e-835a-11eb-364b-09fd2e90b6d8
# ╟─c2f91cf4-835d-11eb-3b5a-fb7f329f4697
# ╠═d4482e58-835d-11eb-1cf1-b1a6d9c5c0b2
# ╟─e1e354b0-835b-11eb-0682-9da5144f8f99
# ╟─0859b1c8-835e-11eb-3be5-2163e163e365
# ╟─46c6770e-831b-11eb-067e-21334c077e86
# ╟─f8ca7706-831c-11eb-1ac3-77fa1483b53f
# ╟─c5b1fb6a-8327-11eb-2e9a-4fc8d2b3eac2
# ╟─58110390-8329-11eb-1c11-8b24865d91ff
# ╠═1c178a14-831d-11eb-2380-93561ebc9bd8
# ╠═dbb19784-831d-11eb-105e-ed8117a41c6b
# ╟─ed24d814-831d-11eb-1018-1d0c85976c0b
# ╟─260bba26-831e-11eb-14b1-47d4ae08397d
# ╟─4107719c-831e-11eb-1766-cb889deb02de
# ╟─c6a576d2-7f85-11eb-1e11-43906c0c1187
# ╟─531698fe-7f85-11eb-0f2d-bbc361c47c72
# ╟─ee13edce-84ed-11eb-394d-85367178ceec
# ╟─4fb1b54a-8320-11eb-0cd4-b38d5e72b1b1
# ╟─a54ba598-831a-11eb-2dcf-dbc008024c55
# ╟─df023fb6-7f70-11eb-215a-7d834c2e7b56
# ╟─06bfbdd0-7c6a-11eb-033b-4df5aaeb3220
# ╠═fbe88cde-7c69-11eb-2c11-f744889098a8
# ╟─9ba79240-7c6c-11eb-3385-7917a37df361
# ╠═74bec6fa-7c6d-11eb-27d0-c972ad74515c
# ╟─803ff050-7c6d-11eb-1bc3-ff5335e4b361
# ╠═14591b22-7c6e-11eb-008f-27a2ca96c8e9
# ╟─29b3d9bc-7c6e-11eb-2702-53b1dddb63ea
# ╠═2111ff28-7c6e-11eb-1ff1-f5017b3889fa
# ╟─955462da-7c6f-11eb-0fe9-a7f3913bf97d
# ╠═bfaf9f4c-7c6f-11eb-01b3-8334e46772b6
# ╟─988d1952-80fb-11eb-3f37-25d71c811413
# ╠═d0a53b18-7c6f-11eb-0afd-0faa5c85a613
# ╟─5f54346a-7866-11eb-3df4-658de2c2ce5f
# ╟─600b116c-84eb-11eb-00e0-f14d7b241d24
# ╟─1d08b996-84eb-11eb-0736-510e5f7cdd1f
# ╟─5cf34e4e-7867-11eb-18ff-01b71b0d5934
# ╟─7c389718-8321-11eb-2e05-d3969090fe0e
# ╟─0926daac-7ccb-11eb-21ca-391422834af3
# ╠═bd2cb4c8-8323-11eb-34dc-a1cdf2cccb25
# ╟─d12e4e0a-8323-11eb-0d74-6b574609a909
# ╠═c54e39f0-7ccb-11eb-3d69-331302b5def2
# ╟─c950647e-7ccb-11eb-028f-d1c7c065342f
# ╠═c89d83a4-7ccb-11eb-210e-355010a04b45
# ╟─e17b1136-7ccb-11eb-1c45-e9de3d9cd831
# ╟─11ea91ce-7cd0-11eb-2aa5-3df9a4b24db5
# ╟─479d8d1c-8324-11eb-2d03-8dfbf74ea283
# ╟─8fdc5676-8324-11eb-01c6-b1f82b890a40
# ╠═ac567efc-7867-11eb-07c3-21989feec15c
# ╟─31795f9a-7b7a-11eb-1786-89aad77aff4b
# ╠═de5fc23a-7b7a-11eb-12b0-a513044b39a6
# ╟─e7475a98-7b7a-11eb-00a2-63e52e403c16
# ╟─e54c99d8-7c74-11eb-29cc-f5d7f64c4937
# ╠═0977c012-7c75-11eb-3cf4-e9a1bfca6541
# ╟─314c5f4e-7c75-11eb-2b48-5787c4c781c9
# ╟─3310e738-7ac1-11eb-38ea-2f1cfc2f0090
# ╠═8f365958-7ac1-11eb-145b-69bff762b7a8
# ╟─9fb90230-7ac1-11eb-1c50-634ba7a9a7bb
# ╠═f3fce35a-7b77-11eb-0935-13399babee72
# ╠═f80d6c4e-7b77-11eb-31bc-d55c9905b0a6
# ╟─fbf35d0a-7b77-11eb-3b36-194f011f00ab
# ╟─ad6d4176-8321-11eb-2cf0-89aa37f483f5
# ╟─0b370ade-7b5d-11eb-14ae-8fb84aa69edc
# ╠═19fd463c-7b5d-11eb-342c-293f27a0f396
# ╠═2f54327a-7b5d-11eb-0352-e713aba2460d
# ╟─3a62156a-7b5d-11eb-1b66-a9d7d71d1051
# ╟─96d3fed8-7b7b-11eb-1434-07da91c27aca
# ╠═c5f3dde8-7b7b-11eb-3a14-bbe7f374c856
# ╟─2a044f08-7b7b-11eb-19fc-9fe29ad589bf
# ╠═8b29e5b4-7b7f-11eb-2679-bd248ed2d88b
# ╟─bc86d990-83d2-11eb-1401-e1ac3ac6b36c
# ╟─95750e6c-83d2-11eb-31e4-673d932bb999
# ╟─0f9373f0-78c5-11eb-26ba-05bea8874141
# ╠═3933c4ee-78c5-11eb-2692-7989afac7ebb
# ╟─6dda76f6-78c6-11eb-0fe0-a36a5663dfd7
# ╠═add65ba2-78c7-11eb-0eb6-771466c0c417
# ╠═d482635e-7ab4-11eb-0d53-398ca84dab54
# ╠═70ab6652-7ab4-11eb-1309-8f833821ef88
# ╟─3cf74998-78c5-11eb-1191-a9f56dbadda6
# ╠═3cfda932-78c5-11eb-3ad9-512459ffebc4
# ╟─3d11b968-78c5-11eb-3962-f50b8e05803d
# ╟─6ff7cb16-83d2-11eb-0155-e38cae990886
# ╟─441ea182-8361-11eb-2fef-11c314784df1
# ╟─fd34d75e-8363-11eb-1ddd-efaf1569a22b
# ╠═e688f190-8361-11eb-1fe8-2bc7c5295926
# ╟─4a4ab3ca-8363-11eb-1c65-cb866fc822ee
# ╠═53b3db4a-8362-11eb-1948-7b49d75d9f99
# ╠═b77ec01c-8363-11eb-037f-833825146fa5
# ╠═75d0456a-8362-11eb-244b-7f53c80e5bc8
# ╠═a7de2550-8363-11eb-038d-696d2e53f139
# ╟─caeae7de-8363-11eb-003d-83d7fbea4217
# ╟─c8daad14-8361-11eb-169d-e36ccbd815ba
# ╟─9fb8ee4a-7cc5-11eb-0a93-2fea93a38c75
# ╟─5974daa4-7cc6-11eb-097a-1303b35124dc
# ╠═c8e07448-7cc6-11eb-1e2d-15d48c74d921
# ╠═d2df623a-7cc6-11eb-2dfd-c57bb3dfddcd
# ╠═df91c5f2-7cc6-11eb-17c4-f9129b2fc701
# ╟─ea286192-7cc6-11eb-3274-6fa2ee43cc44
# ╠═f94f0b80-7cc6-11eb-1ad6-c134a7ae5590
# ╟─03468256-7cc7-11eb-3b5f-a16deac85f21
# ╠═439e7f0e-7cc7-11eb-0080-c970dd4196ba
# ╟─4efe834c-7cc7-11eb-1060-9137f0a8da19
# ╠═64337256-7cc7-11eb-2ef3-4304c425d51d
# ╟─3d198bac-78c5-11eb-35af-87b99956e4de
# ╟─000779f6-7f6c-11eb-012b-6b67b412eb5b
# ╟─1029caaa-7f6c-11eb-0121-5b4b86a9a832
# ╟─5f093552-7f6c-11eb-3ea4-5b0032275023
# ╟─7a843716-7f6c-11eb-3656-a1014b9984ad
# ╟─92dda336-7f6c-11eb-2e80-d959593aff78
# ╟─b8c2920a-7f6c-11eb-37fc-9fa12702d8f2
# ╟─4fb32b18-8330-11eb-1cc8-779e4a09332d
# ╠═63b924f0-8330-11eb-121a-c56a30993c81
# ╟─90a198a8-8330-11eb-3b78-fff35e711256
# ╠═6a25b7e0-8330-11eb-2258-c357b0033d47
# ╠═7a9c88b0-8330-11eb-0938-39c0cca872ae
# ╟─a92a2796-8330-11eb-210b-11fb46621064
# ╠═84e1f6d4-8330-11eb-3fc9-271654090046
# ╟─ef4606b8-8318-11eb-1bf5-5fee8d428bae
# ╟─0698f0c8-8319-11eb-2035-191ff7afc846
# ╟─7bd38294-7b7c-11eb-0329-99e46b366e1b
# ╟─64484bfa-7b7c-11eb-1d40-296fcf167f53
# ╟─0ff8fb0c-7b84-11eb-2a6a-0dedd44b3c7b
# ╟─cc6069f2-7b7c-11eb-3cb9-afbd773708e2
# ╟─18b49d34-7b80-11eb-3b12-3774bdf47ae6
# ╟─b2b0b842-7b57-11eb-1727-95c07b5b2de6
# ╟─9bd103ee-7f75-11eb-0f0d-b13097313f0a
# ╟─b8783bca-7f75-11eb-070c-697f8d79c0e8
# ╟─7779ed4c-7b7d-11eb-034b-5fe2925df5df
