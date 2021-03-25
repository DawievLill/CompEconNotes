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

# ╔═╡ 796b0922-8c17-11eb-31e8-59d5b21ee32b
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
			Pkg.PackageSpec(name="Images"), 
			Pkg.PackageSpec(name="ImageMagick"), 
			Pkg.PackageSpec(name="PlutoUI"), 
			Pkg.PackageSpec(name="HypertextLiteral"), 		
			])

	using Images
	using PlutoUI
	using HypertextLiteral
	using LinearAlgebra
end

# ╔═╡ 97dca378-8c17-11eb-1a9f-49d299180a72
md" # Linear Algebra "

# ╔═╡ 8dd16592-8d9c-11eb-1adb-af30ea9f5bc5
md" In this section we will cover some basic foundational concepts in linear algebra, which will be useful for various types of transformations and eventually solutions of linear systems of equations." 

# ╔═╡ cdcc50b0-8c18-11eb-3188-ffb91dcecee8
md" ## Arrays in Julia "

# ╔═╡ dc8e8004-8d95-11eb-30c0-33a66a0c3f6b
md" An array is a rectangular grid that is used for storing data of any type. The basic array is constructed in a similar way to the that of Matlab, which can be done as follows: "

# ╔═╡ eea539ae-8d95-11eb-1a13-d704f50969c3
x₁ = [1.0, 2.0, 3.0]

# ╔═╡ 1eeb8d5c-8d96-11eb-0bc1-0162c420760f
typeof(x₁) 

# ╔═╡ e9140e2a-8d97-11eb-3a2c-79fad99882ea
md" You might have encountered arrays in mathematics. An array can be one-dimensional, two-dimensional, three-dimensional and so forth. We can see that our example for `x` is a three element array that contains `Float64` values. The `1` means that this is a one-dimensional array, which is also known as a vector. In Julia it is best to think about these as **column vectors**. "

# ╔═╡ 784d31a8-8d9a-11eb-02bc-a3a24bcb564c
typeof(x₁) == Vector{Float64} 

# ╔═╡ ce128ae2-8d9b-11eb-0399-cda8550d4524
md" You will notice if we attempt to create a row vector then it is considered two-dimensional. This means that row vectors are treated as matrices in a sense. This is mainly done for the purposes of multiplication and addition so that dimensions align in the right manner. Multiplying with a column vector is inherently different from multiplying with a row vectors. It is good to take note of these subtleties. " 

# ╔═╡ 98dbd9c8-8d9b-11eb-009d-b93e1b869285
x₂ = [1.0 2.0 3.0]

# ╔═╡ ff747db4-8d9b-11eb-3a10-3314a843e679
typeof(x₂) == Matrix{Float64}

# ╔═╡ b386f1dc-8d9a-11eb-2bb8-a58a56b5a85b
md" We have established that one-dimensional arrays are a vector type. Similarily, we can view a two-dimensional array as a matrix. Consider the following example."

# ╔═╡ cb7dce80-8d9a-11eb-1f83-fb79238361e1
x₃ = [1 2 3; 4 5 6; 7 8 14]

# ╔═╡ 54db62f4-8d9b-11eb-34d6-9f4d16b251de
typeof(x₃) == Matrix{Int64}

# ╔═╡ 1b9b6594-8d9d-11eb-375a-79652da9dee1
ndims(x₃) # Number of dimensions of x₃

# ╔═╡ 2ba16772-8d9d-11eb-0a26-c556fb0240aa
size(x₃) # Tuple containing the dimensions of x₃

# ╔═╡ 859cee86-8d9d-11eb-2993-bb3f3d748df6
length(x₃) # Number of elements in x₃

# ╔═╡ dfa715e2-8d9c-11eb-29e3-7dad0fad127f
md" We see that our $3 \times 3$ matrix is a two-dimensional array. Take note of the difference between dimension and size of the matrix."

# ╔═╡ cd94ab08-8d9c-11eb-2a76-673456a371ef
md" ### Construction of arrays"

# ╔═╡ 9bd1329e-8d9c-11eb-1751-fd348630580e
md" There are many important ways to create and initialise arrays automatically in Julia. Some of the most important methods are shown below. See if you can figure out from the output what the functions are doing.  "

# ╔═╡ ae7906d8-8d9c-11eb-3e39-4fcddb18b1d2
zeros(3)

# ╔═╡ b28764ae-8d9c-11eb-10f7-8f4b2387665f
zeros(3, 3)

# ╔═╡ b7a8d21a-8d9c-11eb-106d-2f68dc163529
fill(5.0, 3, 3)

# ╔═╡ f29ef44a-8d9d-11eb-0dbe-d5cc7c1080db
ones(3, 3)

# ╔═╡ fb532bfe-8d9d-11eb-2a10-89d840f835df
reshape(x₃, 1, 9)

# ╔═╡ 16efc20a-8d9e-11eb-1819-47dd80edfc9b
similar(x₃)

# ╔═╡ 08ef44e4-8da0-11eb-3d23-052a692e711f
rand(Int64, 3, 3)

# ╔═╡ 40b31e3c-8d9e-11eb-0c87-5bf71827d188
rand(x₃, 3, 3) # This one is actually a bit tricky, you need to think about this. 

# ╔═╡ 58232b72-8d9e-11eb-0489-cfccbceb1e11
randn(Float64, 3, 3)

# ╔═╡ dab8d7f8-8d9e-11eb-0805-c7343db6c32b
Matrix{Float64}(I, 3, 3)

# ╔═╡ e5bee5de-8d9e-11eb-12af-617d270d651a
range(0, 100, length = 51) # Similar to `linspace` in Matlab

# ╔═╡ a0d17a54-8c18-11eb-0c42-c1553dfc28d5
md" ## Transformations "

# ╔═╡ a3cecb9e-8c18-11eb-1ce0-b531dfce4c7f
md" In mathematics you will have dealt with matrices and how to multiply them with each other and also with specific scalar values and vectors. Normally we teach you about matrices as tables of numbers while vectors contain columns or rows of values. You then proceed with learning all the rules of multiplication in addition. The good news is that this creates a way of thinking about operations between these mathematical constructs. However, it turns out that we are really bad at calculating these things in practice and that computers can do a much better job. Unfortunately we spend so much time teaching you about the operations of multiplication and addition that we fail to mention the intuition behind these transformations. In this next section we will delve a bit deeper into the intuition behind transformations and let the computer do all the computing. "

# ╔═╡ aff36556-8c18-11eb-1dc1-5d6e8f7f9854
md" #### Transformation of images "

# ╔═╡ b288129e-8c18-11eb-18a8-634cc4c8bba0
img_sources = ["https://user-images.githubusercontent.com/6933510/108883855-39690f80-7606-11eb-8eb1-e595c6c8d829.png" => "Arrows"]

# ╔═╡ b4bdca54-8c18-11eb-081e-2bc3d2ffe410
md"""
You can choose an image from the dropwdown list below. You first need to add the image to the list above, with the appropriate name.

**Choose image from list below:**

$(@bind img_source Select(img_sources))
"""

# ╔═╡ b73be1e4-8c18-11eb-08ab-5bb385271ef5
md" The slider below represents different values of $\alpha$"

# ╔═╡ b8ef5872-8c18-11eb-1683-c32d2aacc4e5
md"""
$(@bind α Slider(.1:.1:3, show_value=true))
"""

# ╔═╡ bad14aa4-8c18-11eb-05f7-55ee2d6e3ad8
let

range = -1.5:.1:1.5
md"""
This is a "scrubbable matrix" -- click on the number and drag to change.	
	
``(``	
 $(@bind a Scrubbable( range; default=1.0))
 $(@bind b Scrubbable( range; default=0.0))
``)``

``(``
$(@bind c Scrubbable(range; default=0.0 ))
$(@bind d Scrubbable(range; default=1.0))
``)``
	
	**Re-run this cell to reset to identity transformation**
"""
end

# ╔═╡ bcfe1640-8c18-11eb-151d-57fb02f34d92
det_A = 1.0

# ╔═╡ bf5de8a4-8c18-11eb-2726-033522175491
begin
	 idy((x,y)) = [x,y]
	 lin1((x,y)) =  [ 2x + 3y, -5x+4x ]
	 scalex(α) = ((x,y),) -> (α*x, y)
	 scaley(α) = ((x,y),) -> (x,   α*y)
	 rot(θ) = ((x,y),) -> [cos(θ)*x + sin(θ)*y, -sin(θ)*x + cos(θ)*y]
	 shear(α) = ((x,y),) -> [x+α*y,y]
	 genlin(a,b,c,d) = ((x,y),) -> [ a*x + b*y ; c*x + d*y ]
end

# ╔═╡ c2fa7130-8c18-11eb-32d3-8b4888e8f7f1
begin
	white(c::RGB) = RGB(1,1,1)
	white(c::RGBA) = RGBA(1,1,1,0.75)
end

# ╔═╡ c5746510-8c18-11eb-251a-37ec2665bf8f
md"""
center zoom = $(@bind z Slider(.1:.1:3, show_value=true, default=1))
"""

# ╔═╡ c14837be-8c18-11eb-207e-2f751082ae79
function trygetpixel(img::AbstractMatrix, x::Float64, y::Float64)
	rows, cols = size(img)
	
	"The linear map [-1,1] ↦ [0,1]"
	f = t -> (t - -1.0)/(1.0 - -1.0)
	
	i = floor(Int, rows *  f(-y) / z)
	j = floor(Int, cols *  f(x * (rows / cols))  / z)
 
	if 1 < i ≤ rows && 1 < j ≤ cols
		img[i,j]
	else
		white(img[1,1])

	end
end

# ╔═╡ c7356bf4-8c18-11eb-0c16-1b5d1bf0fe9b
md"""
top left zoom =	$(@bind f Slider(.1:1:3, show_value=true, default=1))
"""

# ╔═╡ c9f20fa2-8c18-11eb-17e0-4386d34e7387
T = shear(1) # Pick a transformation
#T = genlin(a,b,c,d)

# ╔═╡ cc076990-8c18-11eb-0bd0-a357f823e14d
[
	if det_A == 0
		RGB(1.0, 1.0, 1.0)
	else
		
		 # in_x, in_y = A \ [out_x, out_y]
         # in_x, in_y = xy( [out_x, out_y] )
		in_x, in_y =  T([out_x, out_y])
		trygetpixel(img, in_x, in_y)
	end
	
	for out_y in LinRange(f, -f, 500),
		out_x in LinRange(-f, f, 500)
]

# ╔═╡ Cell order:
# ╟─796b0922-8c17-11eb-31e8-59d5b21ee32b
# ╟─97dca378-8c17-11eb-1a9f-49d299180a72
# ╟─8dd16592-8d9c-11eb-1adb-af30ea9f5bc5
# ╟─cdcc50b0-8c18-11eb-3188-ffb91dcecee8
# ╟─dc8e8004-8d95-11eb-30c0-33a66a0c3f6b
# ╠═eea539ae-8d95-11eb-1a13-d704f50969c3
# ╠═1eeb8d5c-8d96-11eb-0bc1-0162c420760f
# ╟─e9140e2a-8d97-11eb-3a2c-79fad99882ea
# ╠═784d31a8-8d9a-11eb-02bc-a3a24bcb564c
# ╟─ce128ae2-8d9b-11eb-0399-cda8550d4524
# ╠═98dbd9c8-8d9b-11eb-009d-b93e1b869285
# ╠═ff747db4-8d9b-11eb-3a10-3314a843e679
# ╟─b386f1dc-8d9a-11eb-2bb8-a58a56b5a85b
# ╠═cb7dce80-8d9a-11eb-1f83-fb79238361e1
# ╠═54db62f4-8d9b-11eb-34d6-9f4d16b251de
# ╠═1b9b6594-8d9d-11eb-375a-79652da9dee1
# ╠═2ba16772-8d9d-11eb-0a26-c556fb0240aa
# ╠═859cee86-8d9d-11eb-2993-bb3f3d748df6
# ╟─dfa715e2-8d9c-11eb-29e3-7dad0fad127f
# ╟─cd94ab08-8d9c-11eb-2a76-673456a371ef
# ╟─9bd1329e-8d9c-11eb-1751-fd348630580e
# ╠═ae7906d8-8d9c-11eb-3e39-4fcddb18b1d2
# ╠═b28764ae-8d9c-11eb-10f7-8f4b2387665f
# ╠═b7a8d21a-8d9c-11eb-106d-2f68dc163529
# ╠═f29ef44a-8d9d-11eb-0dbe-d5cc7c1080db
# ╠═fb532bfe-8d9d-11eb-2a10-89d840f835df
# ╠═16efc20a-8d9e-11eb-1819-47dd80edfc9b
# ╠═08ef44e4-8da0-11eb-3d23-052a692e711f
# ╠═40b31e3c-8d9e-11eb-0c87-5bf71827d188
# ╠═58232b72-8d9e-11eb-0489-cfccbceb1e11
# ╠═dab8d7f8-8d9e-11eb-0805-c7343db6c32b
# ╠═e5bee5de-8d9e-11eb-12af-617d270d651a
# ╟─a0d17a54-8c18-11eb-0c42-c1553dfc28d5
# ╟─a3cecb9e-8c18-11eb-1ce0-b531dfce4c7f
# ╟─aff36556-8c18-11eb-1dc1-5d6e8f7f9854
# ╟─b288129e-8c18-11eb-18a8-634cc4c8bba0
# ╟─b4bdca54-8c18-11eb-081e-2bc3d2ffe410
# ╟─b73be1e4-8c18-11eb-08ab-5bb385271ef5
# ╟─b8ef5872-8c18-11eb-1683-c32d2aacc4e5
# ╟─bad14aa4-8c18-11eb-05f7-55ee2d6e3ad8
# ╠═bcfe1640-8c18-11eb-151d-57fb02f34d92
# ╟─bf5de8a4-8c18-11eb-2726-033522175491
# ╟─c14837be-8c18-11eb-207e-2f751082ae79
# ╟─c2fa7130-8c18-11eb-32d3-8b4888e8f7f1
# ╠═c5746510-8c18-11eb-251a-37ec2665bf8f
# ╟─c7356bf4-8c18-11eb-0c16-1b5d1bf0fe9b
# ╠═c9f20fa2-8c18-11eb-17e0-4386d34e7387
# ╟─cc076990-8c18-11eb-0bd0-a357f823e14d
