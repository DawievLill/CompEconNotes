

# Bracketing -- Wheeler - p 36 

using Parameters

@with_kw struct params
    s::Float64 = 1e-2
    k::Float64 = 2.0
end

p = params()

function bracket_minimum(f, x)
    @unpack s, k = p

    a, ya = x, f(x)
    b, yb = a + s, f(a + s)

    if yb > ya
        a, b = b, a
        ya, yb = yb, ya
        s = -s
    end 

    while true
        x, yc = b + s, f(b + s)
        if yc > yb
            return a < c ? (a, c) : (c, a)
        end
        a, ya, b, yb = b, yb, c, yc
        s *= k
    end
end
