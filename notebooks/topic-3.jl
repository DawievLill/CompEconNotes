### A Pluto.jl notebook ###
# v0.14.1

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
			Pkg.PackageSpec(name="Random"), 
			Pkg.PackageSpec(name="BenchmarkTools")
			])
	using Random
	using BenchmarkTools
	using Images
	using PlutoUI
	using HypertextLiteral
	using LinearAlgebra
end

# ╔═╡ 97dca378-8c17-11eb-1a9f-49d299180a72
md" # Linear Algebra"

# ╔═╡ 8dd16592-8d9c-11eb-1adb-af30ea9f5bc5
md" In this section we will cover some foundational concepts in (numerical) linear algebra. Our first section focuses on arrays in Julia. Following that we will move on to solving systems of linear equations.  " 

# ╔═╡ d3a32206-96ff-11eb-0c21-953e0d446b50
	md"""
	!!! info
	    For those of you who need a refresher on linear algebra, you can look [here](https://fncbook.github.io/fnc/appendix/linear-algebra.html) or [here](https://julia.quantecon.org/tools_and_techniques/linear_algebra.html). There is also a free book on applied linear algebra that has accompanying Julia code, which you can find [here](http://vmls-book.stanford.edu/). 
	"""

# ╔═╡ cdcc50b0-8c18-11eb-3188-ffb91dcecee8
md" ## Arrays in Julia "

# ╔═╡ dc8e8004-8d95-11eb-30c0-33a66a0c3f6b
md" In this section I will be following similar material to that of the [Quantecon](https://julia.quantecon.org/getting_started_julia/fundamental_types.html) website. In some places the descriptions will be close to (or identical) to their depiction. The purpose of this notebook is to condense all the information into one place. These notes are no substitute for complete coverage on the topic, just a basic introduction.

The first thing we need to consider in this section is our definition for an array. A technical definition is that an array is a rectangular grid that is used for storing data of any type. However, all of the data in the container should be of the same type!

The basic array is constructed in a similar way to the that of Matlab, which can be done as follows: "

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

# ╔═╡ 2f5d84f6-9700-11eb-3dca-799c96c6b609
typeof(size(x₃)) # Immutable type!!

# ╔═╡ 438e6cba-9700-11eb-2c73-339e72e703cb
md" It is important to note here that the *tuple* type is immutable. This means that you cannot change the values once you have created the tuple. Observe the error that we get if we want to alter a specific element of a tuple." 

# ╔═╡ 8f86a65a-9700-11eb-1633-95b7ef7e7295
size(x₃)[1] = 2 # In this case we want to change the first element of the tuple

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
md" The next step is to think about indexing of our constructed matrices / arrays. In other words, how do we get specific values out of our created arrays.  "

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

# ╔═╡ ad4bc944-96f7-11eb-3c46-1f42cb700c68
md" #### Special matrices "

# ╔═╡ 70d37008-9705-11eb-17ac-15f8a398eb86
md" There are some special arrays in Julia that we need take care in using. As always, we need to be careful about types in Julia. Consider the creation of a diagonal matrix." 

# ╔═╡ 9ae75b2a-9705-11eb-1e0e-5d2ca6c0407f
d₁ = [1.0, 2.0]

# ╔═╡ a1cfecfe-9705-11eb-2ea5-2bcdda8c123d
d₂ = Diagonal(d₁) 

# ╔═╡ b62e71e8-9705-11eb-27d1-2b8e7fd1bf4f
md" We observe that d₂ is not a 2-dimensional array as we would expect. The type is something completely different. Julia does this because this type is more efficient to store and to use in arithmetic and matrix operations. We can perform the same operations as before, but be careful to note that this is not a **dense** matrix. "

# ╔═╡ 1146394e-9706-11eb-08be-cd9360651d26
md" Another interesting example is that of the identity matrix. Let us say that we wanted to subtract the identity marix from a 2-dimensional array. We could do it in one of two ways. " 

# ╔═╡ 1e29435e-9706-11eb-1847-9b045e41ac71
d₃ = [1.0 2.0; 3.0 4.0]

# ╔═╡ 3969bd2e-9706-11eb-3c87-8b4ae3261186
d₃ - Diagonal([1.0, 1.0])

# ╔═╡ 7b8586a2-9706-11eb-3062-df72a0a683c9
d₃ - I # Note we didn't specify dimension of I

# ╔═╡ 80e5300c-9706-11eb-393a-6d4dba078567
typeof(I)

# ╔═╡ 8e1b1ffc-9706-11eb-14d3-73f5ae17f967
md" The identity matrix is a strange type that we have not encountered before. It is much more general than the method we employed and also more powerful as it can be applied in many different situations."

# ╔═╡ e846ce68-96f7-11eb-37a5-b706a6979c63
md" ### Basic operations "

# ╔═╡ 4adb0436-9702-11eb-3583-397161d0d7e7
md" Before we talk about systems of linear equations it might be useful to specify some of the basic linear algebra operations that we might encounter."

# ╔═╡ 606ed3a4-9702-11eb-1396-6df331c01343
x₄ = randn(5)

# ╔═╡ 24681522-9703-11eb-33ac-15006eddeb11
y₄ = randn(5)

# ╔═╡ 6a6690ea-9702-11eb-0ee4-0b7ffb0b982b
norm(x₄) # Provides the L2 norm

# ╔═╡ 7f529b0c-9702-11eb-083a-71a10e65776c
sqrt(sum(abs2, x₄)) # This is the L2 norm

# ╔═╡ 96e43f00-9702-11eb-1244-df4c077349c7
md" We might want to also think about calculating dot products. This can be done in the following way: "

# ╔═╡ 51fc430a-9703-11eb-0595-19e9e6a75785
dot(x₄, y₄) # x₄' * y₄

# ╔═╡ ffc56698-9702-11eb-3ab8-0baffe0aa914
x₄'y₄ # Same answer as above

# ╔═╡ b068662e-9703-11eb-13e3-c5898705e274
md" Below is some other useful information that we might want to gather from a matrix. "

# ╔═╡ 7bdfe3e8-9703-11eb-068d-398708020f10
tr(x₃) # Computes the trace of the matrix

# ╔═╡ 97061304-9703-11eb-07f0-b14fbf0fc3e6
det(x₃) # Computes the determinant of the matrix

# ╔═╡ a68db73c-9703-11eb-2bb4-7bbe725c4e3b
rank(x₃) # Computes the rank of the matrix	

# ╔═╡ 6f0c0bc4-96d1-11eb-1cd9-9172bb9f042c
md" # Systems of linear equations"

# ╔═╡ 83ae1458-96d3-11eb-2e21-a3be23f7b874
md" One of the most basic tasks in numerical analysis is to solve the system of linear equations, which can be represented with the following matrix-vector notation. 

$$\textbf{Ax = b}$$

where $\textbf{A}$ is a $n \times n$ matrix, $\textbf{b}$ is a vector of length $n$ and $\textbf{x}$ is a vector of length $n$."

# ╔═╡ 45b2161a-96db-11eb-046a-079f3ddfcee9
md" Judd (1998) states that there are three reasons to solve linear equations.

1. Important class of problems in themselves. 
2. Linear equations are almost the only problem we know how to solve directly. 
3. Many ideas used in this section are applicable to more general problems in later sessions. "

# ╔═╡ 2fbbf6ce-96dc-11eb-28fd-a531a03c2a72
md" For this session we will first discuss **direct methods** for solving linear equations. We follow that with a quick discussion on the condition number for a linear system and then we cover **iterative methods** for solving systems of equations. "

# ╔═╡ 67607a1c-96e3-11eb-3acb-55a989e7efb4
md" ### Direct methods"

# ╔═╡ 2f46eecc-96f5-11eb-3f35-19de6bd828be
md" There are several direct methods for solving linear equations. We will focus on Gaussian elimation (GE) and lower-upper (LU) decomposition. However, one can also look at QR decomposition and singular value decomposition (SVD). I don't believe we will have time for the last two methods, but will perhaps in passing mention how they operate. 

These direct methods are best applied to *dense* and relatively small $\textbf{A}$ matrices."  

# ╔═╡ 86606d6e-96f5-11eb-0991-a7ab19525fde
md" #### Triangular systems "

# ╔═╡ e8ce31d0-9056-48c3-a3d8-0aae4b99eb33
md" We want to solve $\textbf{Ax = b}$. In most cases this means solving $\mathbf{x = A^{-1}b}$. This means that we need to find the inverse of a matrix, which is often hard to do. In the numerical world we **almost never compute an inverse**. The process of finding the solution through an inverse is slower than solving the linear system of equations.

The idea in this section is to the turn the original problem into one that is easy to solve. We would like to transform our $\textbf{A}$ matrix into a structure that can be easily solved. One such representation is the triangular system. Other matrices that are easy to solve include permutation and orthogonal matrices. Below is a representation of a lower triangular matrix.  

$$\left(\begin{array}{cccc}a_{11} & 0 & \cdots & 0 \\ a_{21} & a_{22} & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ a_{n 1} & a_{n 2} & \cdots & a_{n n}\end{array}\right)\left(\begin{array}{c}x_{1} \\ x_{2} \\ \vdots \\ x_{n}\end{array}\right)=\left(\begin{array}{c}b_{1} \\ b_{2} \\ \vdots \\ b_{n}\end{array}\right)$$

The reason we want this strucuture, is that we can then solve by **forward substitution**, which is represented as follows: 

$\begin{aligned} x_{1} &=b_{1} / a_{11} \\ x_{2} &=\left(b_{2}-a_{21} x_{1}\right) / a_{22} \\ x_{3} &=\left(b_{3}-a_{31} x_{1}-a_{32} x_{2}\right) / a_{33} \\ & \quad \vdots \\ x_{n} &=\left(b_{n}-a_{n 1} x_{1}-a_{n 2} x_{2}-\cdots-a_{n, n-1} x_{n-1}\right) / a_{n n} \end{aligned}$"

# ╔═╡ 34ad7ce4-0e9f-40eb-9424-c53ee764e5b9
md" One can code this forward substitution by hand using a `for loop` as follows."

# ╔═╡ f3a02044-f490-4273-8019-8938beba724d
function forwardsub(a,b)

	n = size(a,1)
	x = zeros(n)
	x[1] = b[1]/a[1,1]
	for i = 2:n
    	s = sum( a[i,j]*x[j] for j=1:i-1 )
    	x[i] = ( b[i] - s ) / a[i,i]
	end

	return x
end

# ╔═╡ 120cbbfc-94f3-4187-b8ca-d553e943f9a0
md" Now let us use our function `forwardsub` to solve for a triangular system."

# ╔═╡ fb1c0b79-fd47-4beb-b2cb-0fc671f08164
begin
	Random.seed!(123); # Set seed for reproducibility
	n₁  = 5;
	A₁ = randn(n₁, n₁);
	b₁ = randn(n₁);
end

# ╔═╡ 19af2999-2502-4782-9396-c7d29b2f59cf
A₁

# ╔═╡ 376f2364-7483-4739-8005-bbaac18d649b
L₁ = tril(A₁) # Full triangular matrix

# ╔═╡ 7f85de41-ec85-486f-827f-ea08929e9ab5
md" We can use the `LowerTriangular` command to create a copy of the lower triangular portion of the original $\textbf{A}$ matrix."

# ╔═╡ a15b1ed6-7270-4feb-8fcc-498249e5ff91
L₂ = LowerTriangular(A₁)

# ╔═╡ c0290bf7-96c7-42ef-80d1-35f6f88f08eb
b₁

# ╔═╡ 6a23242c-ed7f-4f45-a1fa-3289f86f7156
md" We have written our function to solve this problem. Now given that you have a triangular structure, how would you go about solving this with some of the included functions in Julia?  Using the `\` operator in Julia will solve the problem using forward substitution. It detects that there is a lower triangular structure and dispacthes to a basic linear algebra subprogram (BLAS) routine to compute. In this case the subroutine is most likely [trsv](https://www.netlib.org/lapack/explore-html/d6/d96/dtrsv_8f.html). Let us compare the two methods."

# ╔═╡ ba9515e2-eb87-441c-ad54-ab018fc74625
forwardsub(L₂, b₁)

# ╔═╡ c5188bbe-1d2f-4964-aa99-0a151118c4af
L₂ \ b₁

# ╔═╡ ff3d034c-fd35-47e1-819a-c048ea5f2ec4
md" We get the same answer. So is the Julia routine faster?"

# ╔═╡ b8fb29ee-f244-4060-ba85-6d0518426691
@benchmark forwardsub(L₂, b₁)

# ╔═╡ e279a91e-5fab-4b73-8f63-610750c2fd1d
@benchmark L₂ \ b₁

# ╔═╡ af09ca7e-fcab-48fe-b61b-115edaa814d7
md" #### Gaussian elimination " 

# ╔═╡ 68582594-c03b-4c1b-ad53-3549cedc1e8b
md" We have encountered the `\` operator before in the case of a triangular system. However, what happens when we call this operator on solve a linear equation when the system is not triangular? In other words, what happens when we call `A \ b` to solve a linear equation? "

# ╔═╡ 6b82fef2-e19a-4551-88b4-6093d6d7ff7f
A₂ = [2.0 1.0 -1.0; -3.0 -1.0 2.0; -2.0 1.0 2.0]

# ╔═╡ f29371ee-4e2e-48bc-95a6-5c52820504ba
b₂ = [8.0, -11.0, -3.0]

# ╔═╡ b0dd755c-3793-44fe-b6d9-fce3b3289564
A₂ \ b₂ 

# ╔═╡ f89e0088-dbc4-43df-8159-39184b735930
md" In order to understand the process of Gaussian elimination, you have to know what an elementary matrix represents. For those who know the procedure feel free to move on to the example. An elementary matrix is an identity matrix with the zero in position $(j, k)$ replaced by a value of $c$. 

$\mathbf{E}_{j k}(c)=\left(\begin{array}{ccccccc}1 & & & & & & \\ & \ddots & & & & & \\ & & 1 & & & & \\ & & & \ddots & & & \\ & & c & & 1 & & \\ & & & & & \ddots & \\ & & & & & & 1\end{array}\right)=\mathbf{I}+c \mathbf{e}_{j} \mathbf{e}_{k}^{T}$


The term $\mathbf{E}_{j k}(c)$ left multiplies a matrix $\mathbf{X}$ to replace its $j$-th row $\mathbf{x}_{j}$ by $c\mathbf{x_k} + \mathbf{x_j}$:

$\mathbf{E}_{j k}(c) \times \mathbf{X}=\mathbf{E}_{j k}(c) \times\left(\begin{array}{llll}\ldots & \mathbf{x}_{k} & \ldots \\ \cdots & \mathbf{x}_{j} & \cdots\end{array}\right)=\left(\begin{array}{lll}\ldots & \mathbf{x}_{k \cdot} & \ldots \\ \cdots & c \mathbf{x}_{k}+\mathbf{x}_{j} & \ldots\end{array}\right)$

Gaussian elimination effectively applies a sequence of elementary matrices to transform the linear system that we are provided to an upper triangular structure $\mathbf{U}$.

$$\begin{aligned} \mathbf{E}_{n, n-1}\left(c_{n, n-1}\right) \cdots \mathbf{E}_{21}\left(c_{21}\right) \mathbf{A} \mathbf{x} &=\mathbf{E}_{n, n-1}\left(c_{n, n-1}\right) \cdots \mathbf{E}_{21}\left(c_{21}\right) \mathbf{b} \\ \mathbf{U x} &=\mathbf{b}_{\mathrm{new}} \end{aligned}$$"

# ╔═╡ 10e36444-e747-4cde-a51c-1e57d3a6cb9e
md" To illustrate the procedure, let us show the steps in the proces. For the first column we construct our elementary matrix with the hope of eliminating the `-3.0` in the `(2, 1)` position of the $A_2$ matrix. Using a value of `1.5` in the `(2, 1)` position of the $E_{21}$ elementary matrix provides us with the required result, as we show below. Can you figure out why we use a value of `1.5`?"

# ╔═╡ 6e278a30-363d-413b-8c7e-833b9b0d1047
E21 = [1.0 0.0 0.0; 1.5 1.0 0.0; 0.0 0.0 1.0]

# ╔═╡ d328e1eb-4feb-4e40-9835-bf5e64beb263
md" We left multiply $E_{21}$ with $A_2$ to replace its second row `[-3.0, -1.0, 2.0]` by `1.5*[2.0, 1.0, -1.0] + [-3.0, -1.0, 2.0]`. This will give us the following."

# ╔═╡ fba688bb-a451-450b-bc1a-9002e090412c
E21 * A₂ # This gives us a zero in the (2, 1) position.

# ╔═╡ 30216472-4160-459d-8796-0b2fa1e9b20e
md" We still need a zero in the `(3, 1)` and `(3, 2)` position. All steps for the process are shown below. " 

# ╔═╡ d73e2a8d-f2df-4353-a4a8-141d34af93e8
E31 = [1.0 0.0 0.0; 0.0 1.0 0.0; 1.0 0.0 1.0]

# ╔═╡ 3396b7ff-2022-4bda-a5cc-98d330f6b103
E31 * E21 * A₂ # This gives us a zero in the (2, 1) and (3, 1) position.

# ╔═╡ 222532b1-fdc3-4aa8-ae40-48fdc54c9e26
E32 = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 -4.0 1.0]

# ╔═╡ b6cc62b8-1f4b-4fdc-9376-3365908cf428
E32 * E31 * E21 * A₂ # We have our upper triangular structure. 

# ╔═╡ 5a3e13a9-78a6-46e0-bec5-5fa18b9db13d
md" We performed our transformation one column at a time. If we combine the elementary matrix operations for each column into one operator, that is referred to as a Gauss transformation. In general for the first column the formula would be: 

$\mathbf{M}_{1}=\mathbf{E}_{n 1}\left(c_{n 1}\right) \cdots \mathbf{E}_{31}\left(c_{31}\right) \mathbf{E}_{21}\left(c_{21}\right)$

In our case we have that $\mathbf{M}_{1} = \mathbf{E}_{31}\left(c_{31}\right) \mathbf{E}_{21}\left(c_{21}\right)$ and $\mathbf{M}_{2} = \mathbf{E}_{32}\left(c_{32}\right)$."

# ╔═╡ 153453d8-2de0-4a40-86b9-24cf83e2d8ee
M1 = E31 * E21

# ╔═╡ 8ee0c629-8c1e-4d61-b9fe-d2fe3a9543ed
M2 = E32

# ╔═╡ 3e1db869-94ab-4db0-8daa-0bacd737e408
md" In essence Gaussian elimination does $\mathbf{M}_{n-1} \cdots \mathbf{M}_{1} \mathbf{A}=\mathbf{U}$. Next time we will continue with the LU decomposition, but I think this should be enough to keep you busy for a while. "

# ╔═╡ bbfe4b6e-8b72-4f6b-b516-1b8838adbcee
md" #### LU factorisation "

# ╔═╡ 7573a640-96e3-11eb-1214-070209074966
md" ### Norms and condition numbers "

# ╔═╡ 8139e91c-96e3-11eb-1d43-7d9502ac6d91
md" ### Iterative solvers "

# ╔═╡ ce083a66-96f5-11eb-1e1c-639e4764cc51
md" There are many different iterative methods, we will focus on the Jacobi method, Gauss-Seidel method and successive over-relaxation (SOR) in this section. Conjugate-gradient methods will also be briefly mentioned. 

These iterative methods are best applied to large, *sparse*, structured linear systems." 

# ╔═╡ Cell order:
# ╟─796b0922-8c17-11eb-31e8-59d5b21ee32b
# ╟─97dca378-8c17-11eb-1a9f-49d299180a72
# ╟─8dd16592-8d9c-11eb-1adb-af30ea9f5bc5
# ╟─d3a32206-96ff-11eb-0c21-953e0d446b50
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
# ╠═2f5d84f6-9700-11eb-3dca-799c96c6b609
# ╟─438e6cba-9700-11eb-2c73-339e72e703cb
# ╠═8f86a65a-9700-11eb-1633-95b7ef7e7295
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
# ╟─ad4bc944-96f7-11eb-3c46-1f42cb700c68
# ╟─70d37008-9705-11eb-17ac-15f8a398eb86
# ╠═9ae75b2a-9705-11eb-1e0e-5d2ca6c0407f
# ╠═a1cfecfe-9705-11eb-2ea5-2bcdda8c123d
# ╟─b62e71e8-9705-11eb-27d1-2b8e7fd1bf4f
# ╟─1146394e-9706-11eb-08be-cd9360651d26
# ╠═1e29435e-9706-11eb-1847-9b045e41ac71
# ╠═3969bd2e-9706-11eb-3c87-8b4ae3261186
# ╠═7b8586a2-9706-11eb-3062-df72a0a683c9
# ╠═80e5300c-9706-11eb-393a-6d4dba078567
# ╟─8e1b1ffc-9706-11eb-14d3-73f5ae17f967
# ╟─e846ce68-96f7-11eb-37a5-b706a6979c63
# ╟─4adb0436-9702-11eb-3583-397161d0d7e7
# ╠═606ed3a4-9702-11eb-1396-6df331c01343
# ╠═24681522-9703-11eb-33ac-15006eddeb11
# ╠═6a6690ea-9702-11eb-0ee4-0b7ffb0b982b
# ╠═7f529b0c-9702-11eb-083a-71a10e65776c
# ╟─96e43f00-9702-11eb-1244-df4c077349c7
# ╠═51fc430a-9703-11eb-0595-19e9e6a75785
# ╠═ffc56698-9702-11eb-3ab8-0baffe0aa914
# ╟─b068662e-9703-11eb-13e3-c5898705e274
# ╠═7bdfe3e8-9703-11eb-068d-398708020f10
# ╠═97061304-9703-11eb-07f0-b14fbf0fc3e6
# ╠═a68db73c-9703-11eb-2bb4-7bbe725c4e3b
# ╟─6f0c0bc4-96d1-11eb-1cd9-9172bb9f042c
# ╟─83ae1458-96d3-11eb-2e21-a3be23f7b874
# ╟─45b2161a-96db-11eb-046a-079f3ddfcee9
# ╟─2fbbf6ce-96dc-11eb-28fd-a531a03c2a72
# ╟─67607a1c-96e3-11eb-3acb-55a989e7efb4
# ╟─2f46eecc-96f5-11eb-3f35-19de6bd828be
# ╟─86606d6e-96f5-11eb-0991-a7ab19525fde
# ╟─e8ce31d0-9056-48c3-a3d8-0aae4b99eb33
# ╟─34ad7ce4-0e9f-40eb-9424-c53ee764e5b9
# ╠═f3a02044-f490-4273-8019-8938beba724d
# ╟─120cbbfc-94f3-4187-b8ca-d553e943f9a0
# ╠═fb1c0b79-fd47-4beb-b2cb-0fc671f08164
# ╠═19af2999-2502-4782-9396-c7d29b2f59cf
# ╠═376f2364-7483-4739-8005-bbaac18d649b
# ╟─7f85de41-ec85-486f-827f-ea08929e9ab5
# ╠═a15b1ed6-7270-4feb-8fcc-498249e5ff91
# ╠═c0290bf7-96c7-42ef-80d1-35f6f88f08eb
# ╟─6a23242c-ed7f-4f45-a1fa-3289f86f7156
# ╠═ba9515e2-eb87-441c-ad54-ab018fc74625
# ╠═c5188bbe-1d2f-4964-aa99-0a151118c4af
# ╟─ff3d034c-fd35-47e1-819a-c048ea5f2ec4
# ╠═b8fb29ee-f244-4060-ba85-6d0518426691
# ╠═e279a91e-5fab-4b73-8f63-610750c2fd1d
# ╟─af09ca7e-fcab-48fe-b61b-115edaa814d7
# ╟─68582594-c03b-4c1b-ad53-3549cedc1e8b
# ╠═6b82fef2-e19a-4551-88b4-6093d6d7ff7f
# ╠═f29371ee-4e2e-48bc-95a6-5c52820504ba
# ╠═b0dd755c-3793-44fe-b6d9-fce3b3289564
# ╟─f89e0088-dbc4-43df-8159-39184b735930
# ╟─10e36444-e747-4cde-a51c-1e57d3a6cb9e
# ╠═6e278a30-363d-413b-8c7e-833b9b0d1047
# ╟─d328e1eb-4feb-4e40-9835-bf5e64beb263
# ╠═fba688bb-a451-450b-bc1a-9002e090412c
# ╟─30216472-4160-459d-8796-0b2fa1e9b20e
# ╠═d73e2a8d-f2df-4353-a4a8-141d34af93e8
# ╠═3396b7ff-2022-4bda-a5cc-98d330f6b103
# ╠═222532b1-fdc3-4aa8-ae40-48fdc54c9e26
# ╠═b6cc62b8-1f4b-4fdc-9376-3365908cf428
# ╟─5a3e13a9-78a6-46e0-bec5-5fa18b9db13d
# ╠═153453d8-2de0-4a40-86b9-24cf83e2d8ee
# ╠═8ee0c629-8c1e-4d61-b9fe-d2fe3a9543ed
# ╟─3e1db869-94ab-4db0-8daa-0bacd737e408
# ╟─bbfe4b6e-8b72-4f6b-b516-1b8838adbcee
# ╟─7573a640-96e3-11eb-1214-070209074966
# ╟─8139e91c-96e3-11eb-1d43-7d9502ac6d91
# ╟─ce083a66-96f5-11eb-1e1c-639e4764cc51
