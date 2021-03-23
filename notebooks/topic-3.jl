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

# ╔═╡ f9dfcf7a-8c17-11eb-00e5-f3ff9b502b4b
md" In this section we will cover some basic foundational concepts in linear algebra, which will be useful for various types of transformations and eventually solutions of linear systems of equations." 

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

# ╔═╡ cdcc50b0-8c18-11eb-3188-ffb91dcecee8
md" ### Arrays in Julia "

# ╔═╡ cf5d7276-8c18-11eb-0dfe-bf95183d1f7d
md" An array is a rectangular grid that is used for storing data of any type. We can store and retrieve data within this array by using indexing. You might have encountered arrays in mathematics under the name of matrices. An array can be one-dimensional, two-dimensional, three-dimensional and so forth. The dimension gives us an idea of the number of indices that we need to specify. For array objects we also need to know the length of the data in each of the dimensions. "

# ╔═╡ d0db00f8-8c18-11eb-13a3-6908bc83cede
md" It might sound arbitrary at first to focus on things like dimension and length of arrays, but take it from me, you will frequently want to know the properties of the arrays that you are working with."

# ╔═╡ Cell order:
# ╟─796b0922-8c17-11eb-31e8-59d5b21ee32b
# ╟─97dca378-8c17-11eb-1a9f-49d299180a72
# ╟─f9dfcf7a-8c17-11eb-00e5-f3ff9b502b4b
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
# ╟─cdcc50b0-8c18-11eb-3188-ffb91dcecee8
# ╟─cf5d7276-8c18-11eb-0dfe-bf95183d1f7d
# ╟─d0db00f8-8c18-11eb-13a3-6908bc83cede
