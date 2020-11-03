



function fibonacci_search(f, a, b, n; ϵ=0.01)

    s = (1 - √5)/(1 + √5)
    ρ = 1 / (φ * (1 - s^(n + 1))/(1 - s^n))
    d = ρ * b + (1 - ρ) * a
    yd = f(d)

    for i in 1 : n-1
        if i == n-1
            c = ϵ * a + (1 - ϵ) * d
        else
            c = ρ * a + (1 - ρ) * b
        end
        yc = f(c)
        if yc < yd
            b, d, yd = d, c, yc
        else
            a, b = b, c
        end
        ρ = 1 / (φ * (1 - s^(n - i + 1))/(1 - s^(n - i)))
    end

    return a < b ? (a, b) : (b, a)

end
    