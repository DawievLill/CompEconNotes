
# Newton-Cotes trapezoid

"""
trapezoid(f,a,b,n)

Apply the trapezoid integration formula for integrand `f` over interval [`a`,`b`], broken up into `n` equal pieces. Returns estimate, vector of nodes, and vector of integrand values at the nodes.
"""
function trapezoid(f,a,b,n)
    h = (b-a)/n
    t = LinRange(a,b,n+1)
    y = f.(t)
    T = h * ( sum(y[2:n]) + 0.5*(y[1] + y[n+1]) )

    return T, t, y
end