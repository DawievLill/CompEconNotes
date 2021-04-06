### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

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
md" In this section we will cover some basic foundational concepts in linear algebra, which will be useful for various types of transformations and eventually solutions of linear systems of equations. Our discussion on transformations follows in a sense from our discussion on functions from last time, where a transformation is a mapping from multiple inputs to multiple outputs. " 

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

# ╔═╡ a1cd647e-92f1-11eb-0ddb-fdc6a622bb8e
size(x₃, 1) # Number of rows

# ╔═╡ af513998-92f1-11eb-16c2-3bd34fe8a231
size(x₃, 2) # Number of columns

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

# ╔═╡ e6c21a46-9311-11eb-3655-0dce4406e2c5
md" The next logical step is to think about indexing of our constructed matrices / arrays. "

# ╔═╡ f03aece2-9311-11eb-1467-c73576cd3184
md" ### Indexing "

# ╔═╡ 1398128c-9312-11eb-3f28-07be67646964
md" The indexing will be very familiar to those that use Matlab and R, not so much for Python. Let us work with x₃ here. Recall that that is the following matrix. "

# ╔═╡ ff94ad36-9311-11eb-297d-5192bd1882ec
x₃

# ╔═╡ 113271b0-9312-11eb-13c2-cd42aaf2ba53
x₃[:, 1] # Provides us with the first column (all the rows of the first column) 

# ╔═╡ 60e1fce2-9312-11eb-19db-df30cb82005f
x₃[1, :] # Provides the first row (all the columns of the first row)

# ╔═╡ 44e3c3c2-9312-11eb-34a4-e3a6b4e5d937
x₃[1:2, 2:3] # Creates a sub-array. Can you figure out how this works?

# ╔═╡ f2d1cce0-9312-11eb-0747-3f0c6f75c1c2
md" **Important!** Creating subsets of an array creates a copy in Julia! It does not alter the original array. See what happens in the following example: " 

# ╔═╡ 35d91f66-9313-11eb-075c-3ba0a8b38f60
z₃ = x₃[1:2, 2:3]

# ╔═╡ 47bd69a6-9313-11eb-0aff-0116148e2a69
z₃[2, 2] = 1.0

# ╔═╡ 551fc65e-9313-11eb-3043-516e2886485f
z₃

# ╔═╡ 59639cac-9313-11eb-008e-9fd2113271e7
x₃ # x₃ did not change as a result of the change in z₃. 

# ╔═╡ 6fb12d46-9313-11eb-3755-f75318201c67
md" This is not standard across programming languages. In other languages the operations performed will alter the original x₃. In order to do this in Julia, we have to use the `views` macro." 

# ╔═╡ 98ada18e-9313-11eb-1b00-a7df44a94629
@views y₃ = x₃[1:2, 2:3]

# ╔═╡ b4b7fe88-9313-11eb-1360-eda276fd3048
y₃[2, 2] = 1.0

# ╔═╡ c2ec5cba-9313-11eb-3a7f-45268de96491
y₃

# ╔═╡ c63186d4-9313-11eb-0c5c-f93db7a91c3d
x₃ # Look at the value of `1` in the 2nd row and third column. This is different from our original matrix. 

# ╔═╡ 44cd0f9a-9314-11eb-375e-2994eee141b8
md" Now for the somewhat confusing part. In the case of subsetting we have that a copy is created, but if you use equality between variables (arrays), then we have that both variables point to the same data. They are allocated to the same place in memory. "

# ╔═╡ 81cfcd24-9314-11eb-32ae-55253c36edee
w₃ = x₃

# ╔═╡ 98f43312-9314-11eb-3442-33b1fecc4542
pointer(w₃), pointer(x₃) # Points to same data

# ╔═╡ b6647ccc-9314-11eb-3c5e-236ba616007e
md" If we were to change w₃ in any way, it would also change x₃." 

# ╔═╡ cd203af2-9314-11eb-08a8-a5eb297100d2
w₃[:, 1] .= 0

# ╔═╡ dd873382-9314-11eb-2890-038d9bb37b5a
x₃

# ╔═╡ e684a78a-9314-11eb-01cf-d3f8eb0b9129
md" If we want to create a copy of x₃ so that we can work with the newly created variable and not alter the original value of x₃, then we can use the `copy` function."

# ╔═╡ 12a4536a-9315-11eb-0c2e-0dbf2ba21af1
v₃ = copy(x₃)

# ╔═╡ 233e7d4a-9315-11eb-38a0-19744474c3e7
pointer(v₃), pointer(x₃) # Pointing to different data

# ╔═╡ 37862b4c-9315-11eb-3e9e-b36398d242ee
md" The general lesson is that subsetting creates a `copy`, and setting arrays equal to each other creates a `view`. If you are aware of these then you won't have a problem. "

# ╔═╡ 6f0c0bc4-96d1-11eb-1cd9-9172bb9f042c
md" # Systems of linear equations"

# ╔═╡ 83ae1458-96d3-11eb-2e21-a3be23f7b874
md" One of the most basic tasks in numerical analysis is to solve the system of linear equations, which can be represented with the following matrix-vector notation. 

$$\textbf{Ax = b}$$

where $\textbf{A}$ is a $m \times n$ matrix, $\textbf{b}$ is a vector of length $m$ and $\textbf{x}$ is a vector of length $n$."

# ╔═╡ 45b2161a-96db-11eb-046a-079f3ddfcee9
md" Judd (1998) states that there are three reasons to solve linear equations.

1. Important class of problems in themselves. 
2. Linear equations are almost the only problem we know how to solve directly. 
3. Many ideas used in this section are applicable to more general problems in later sessions. "

# ╔═╡ 2fbbf6ce-96dc-11eb-28fd-a531a03c2a72
md" For this session we will first discuss **direct methods** for solving linear equations. We quickly discuss the issue of a condition number for a linear system and then we cover **iterative methods** for solving systems of equations. "

# ╔═╡ a0d17a54-8c18-11eb-0c42-c1553dfc28d5
md" ## Transformations "

# ╔═╡ a3cecb9e-8c18-11eb-1ce0-b531dfce4c7f
md" In mathematics you will have dealt with matrices and how to multiply them with each other and also with specific scalar values and vectors. Normally we teach you about matrices as tables of numbers while vectors contain columns or rows of values. You then proceed with learning all the rules of multiplication in addition. The good news is that this creates a way of thinking about operations between these mathematical constructs. However, it turns out that we are really bad at calculating these things in practice and that computers can do a much better job. Unfortunately we spend so much time teaching you about the operations of multiplication and addition that we fail to mention the intuition behind these transformations. In this next section we will delve a bit deeper into the intuition behind transformations and let the computer do all the computing. "

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
# ╠═a1cd647e-92f1-11eb-0ddb-fdc6a622bb8e
# ╠═af513998-92f1-11eb-16c2-3bd34fe8a231
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
# ╟─e6c21a46-9311-11eb-3655-0dce4406e2c5
# ╟─f03aece2-9311-11eb-1467-c73576cd3184
# ╟─1398128c-9312-11eb-3f28-07be67646964
# ╠═ff94ad36-9311-11eb-297d-5192bd1882ec
# ╠═113271b0-9312-11eb-13c2-cd42aaf2ba53
# ╠═60e1fce2-9312-11eb-19db-df30cb82005f
# ╠═44e3c3c2-9312-11eb-34a4-e3a6b4e5d937
# ╟─f2d1cce0-9312-11eb-0747-3f0c6f75c1c2
# ╠═35d91f66-9313-11eb-075c-3ba0a8b38f60
# ╠═47bd69a6-9313-11eb-0aff-0116148e2a69
# ╠═551fc65e-9313-11eb-3043-516e2886485f
# ╠═59639cac-9313-11eb-008e-9fd2113271e7
# ╟─6fb12d46-9313-11eb-3755-f75318201c67
# ╠═98ada18e-9313-11eb-1b00-a7df44a94629
# ╠═b4b7fe88-9313-11eb-1360-eda276fd3048
# ╠═c2ec5cba-9313-11eb-3a7f-45268de96491
# ╠═c63186d4-9313-11eb-0c5c-f93db7a91c3d
# ╟─44cd0f9a-9314-11eb-375e-2994eee141b8
# ╠═81cfcd24-9314-11eb-32ae-55253c36edee
# ╠═98f43312-9314-11eb-3442-33b1fecc4542
# ╟─b6647ccc-9314-11eb-3c5e-236ba616007e
# ╠═cd203af2-9314-11eb-08a8-a5eb297100d2
# ╠═dd873382-9314-11eb-2890-038d9bb37b5a
# ╟─e684a78a-9314-11eb-01cf-d3f8eb0b9129
# ╠═12a4536a-9315-11eb-0c2e-0dbf2ba21af1
# ╠═233e7d4a-9315-11eb-38a0-19744474c3e7
# ╟─37862b4c-9315-11eb-3e9e-b36398d242ee
# ╟─6f0c0bc4-96d1-11eb-1cd9-9172bb9f042c
# ╟─83ae1458-96d3-11eb-2e21-a3be23f7b874
# ╟─45b2161a-96db-11eb-046a-079f3ddfcee9
# ╟─2fbbf6ce-96dc-11eb-28fd-a531a03c2a72
# ╟─a0d17a54-8c18-11eb-0c42-c1553dfc28d5
# ╟─a3cecb9e-8c18-11eb-1ce0-b531dfce4c7f
