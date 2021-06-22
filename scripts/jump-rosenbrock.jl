

# JuMP: nonlinear Rosenbrock Example
# Instead of hand-coding first and second derivatives, you only have to give `JuMP` expressions for objective and constraints.
# Here is an example.

using JuMP
using Ipopt

let

    m₁ = Model(with_optimizer(Ipopt.Optimizer))

    @variable(m₁, x₁)
    @variable(m₁, y₁)

    @NLobjective(m₁, Min, (1-x₁)^2 + 100(y₁-x₁^2)^2)

    JuMP.optimize!(m₁)
    @show value(x₁)
    @show value(y₁)
    @show termination_status(m₁)

end

# not bad, right?
# adding the constraint from before:

let
    
    m₂ = Model(with_optimizer(Ipopt.Optimizer))

    @variable(m₂, x₂)
    @variable(m₂, y₂)

    @NLobjective(m₂, Min, (1-x₂)^2 + 100(y₂-x₂^2)^2)


    @NLconstraint(m₂, con, x₂^2 + y₂^2 <= 0.8)

    JuMP.optimize!(m₂)
    @show value(x₂)
    @show value(y₂)
    @show termination_status(m₂)
    @show dual(con)

end