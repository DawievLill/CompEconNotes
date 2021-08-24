
# Live coding session

using LinearAlgebra

A = I(2)
cond(A)

ϵ = 1E-6

A2 = [1.0 0.0 
     1.0 ϵ]

cond(A2)     

inv(A2)

det(A2) # poor conditioning

# Gauss-Jacobi, Gauss-Seidel, Conjugate gradient (newest) 


