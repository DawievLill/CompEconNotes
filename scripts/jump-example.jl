
# Simple example using JuMP

using JuMP
using GLPK

model = Model(with_optimizer(GLPK.Optimizer)) # Specify the optimizer
@variable(model, 0 <= x <= 2) # Provide range for x variable 
@variable(model, 0 <= y <= 30) # Provide range for y variable 

# next, we set an objective function
@objective(model, Max, 5x + 3y)

# maybe add a constraint called "con"?
@constraint(model, con, 1x + 5y <= 3)

JuMP.optimize!(model)

# look at status
termination_status(model)

# we query objective value and solutions
@show objective_value(model)
@show value(x)
@show value(y)

# as well as the value of the dual variable on the constraint (Lagrange multiplier)
@show dual(con);
