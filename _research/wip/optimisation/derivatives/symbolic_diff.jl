
# Wheeler - p21 (Derivatives)

# Limit equation defining the derivative can be presented in three different ways, forward difference, central difference and backward difference

# If the function can be represented symbolically we can use symbolic differentiation via the SymEngine package in Julia -- Sympy in python. 

using SymEngine

@vars x; # This defines x as a symbolic variable
f = x^2 + x/2 - sin(x)/x;
diff(f, x) # Symbolic differentiation