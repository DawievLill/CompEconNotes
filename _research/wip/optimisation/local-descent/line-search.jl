#=
Line search from local descent chapter

@author: Dawie van Lill <dvanlill@sun.ac.za>

@date: 2020-11-09

References
----------
Algorithms for Optimisation

=#

function bracket_minimum(f, x=0; s=1e-2, k=2.0)
    a, ya = x, f(x)
    b, yb = a + s, f(a + s)
    
    if yb > ya
        a, b = b, a
        ya, yb = yb, ya
        s = -s
    end
    
    while true
        c, yc = b + s, f(b + s)
        if yc > yb
            return a < c ? (a, c) : (c, a)
        end
        a, ya, b, yb = b, yb, c, yc
        s *= k
    end
end

function bracket_sign_change(f′, a, b; k=2)
    if a > b; a,b = b,a; end # ensure a < b
    center, half_width = (b+a)/2, (b-a)/2
    
    while f′(a)*f′(b) > 0
        half_width *= ka
        a = center - half_width
        b = center + half_width
    end
    
    return (a,b)
end
    

function line_search(f, x, d)
    objective = α -> f(x + α * d)
    a, b = bracket_minimum(objective)
    α = bracket_sign_change(objective, a, b)

    return x + α * d
end