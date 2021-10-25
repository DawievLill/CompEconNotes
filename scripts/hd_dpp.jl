
# Julia replication of Solving High-Dimensional Dynamic Programming Problems using Deep Learning
# Original code written for Tensorflow by Jesús Fernández-Villaverde (University of Pennsylvania, NBER and CEPR), Roberto-Rafael Maura-Rivero (London School of Economics), Galo Nuño (Banco de España), George Sorg-Langhans (Princeton University) and Maximilian Vogler (Princeton University)

# Required packages
using Distributions
using Flux
using LinearAlgebra
using Parameters
using Plots
using Random
using Zygote


# Parameter defaults (MWE, so not all included)
@with_kw struct Params
    
    # Note: type declarations not strictly required (could have impact on speed, but unlikely)

    γ::Float64 = 2.0 # curvature of the utility function
    ρ::Float64 = 0.04 # discount rate
    A::Float64 = 0.5 # total factor productivity
    α::Float64 = 0.36 # returns to scale
    δ::Float64 = 0.05 # depreciation rate of capital

    k_min::Float64 = 0.1 
    k_max::Float64 = 10.0
    batch_size::Int64 = 1000

    # ergodic distribution estimation (not used in this setup)
    n_burned::Int64 = 1 # number of periods iterated to move away from starting point

    # parameters for neural network
    n_neurons::Int64 = 8
    input_dim::Int64 = 1

    # misc
    dt::Float64 = 0.1
    epsilon::Int64 = 1
    
end

p = Params() # change parameter defaults here

### MODEL

# Value function neural net
function neural_net_vf(p::Params; bias = true)
    @unpack n_neurons, input_dim = p

    # neural network with three hidden layers
    return model = Chain(
        Dense(input_dim, n_neurons, tanh), 
        Dense(n_neurons, n_neurons, tanh),
        Dense(n_neurons, n_neurons, tanh),
        Dense(n_neurons, 1, bias = -60 * ones(1))
    )
end

# Policy function neural net
function neural_net_cf(p::Params; bias = true)
    @unpack n_neurons, input_dim = p

    # neural network with three hidden layers
    return model = Chain(
        Dense(input_dim, n_neurons, tanh), 
        Dense(n_neurons, n_neurons, tanh),
        Dense(n_neurons, n_neurons, tanh),
        Dense(n_neurons, 1, bias = zeros(1))
    )
end

# This generates 169 parameters, which is the same as the TF code

#### HJB EQUATION

# Value error stems from this
function hjb(p::Params, k_capital, V, C)
    @unpack A, α, δ, ρ, γ = p

    v_prime = Zygote.gradient(k_capital -> sum(V(k_capital)), k_capital)[1] # Have to take the sum to reflect same as TensorFlow

    Y = A * k_capital .^ α # output
    I = Y - exp.(C(k_capital)) # investment
    dK_dt = I - δ * k_capital # capital drift (dK/dt)
    
    # check the explicit functional form from the paper
    U = (exp.(C(k_capital)) .^ (1 - γ)) / (1 - γ) # utility 
    hjb = U - ρ * V(k_capital) + Zygote.dropgrad(v_prime) .* dK_dt 

    return hjb
end

function ergodic_distribution(p::Params, k_capital, C)
    @unpack A, α, δ, ρ, n_burned, dt = p

    # iterate points to estimate the ergodic distribution (check the range of the for loop here)
    for i in 1:n_burned

        Y = A * k_capital .^ α # output
        I = Y - exp.(C(k_capital)) # investment
        dK_dt = I - δ * k_capital # capital drift 

        k_capital = k_capital + (dt .* dK_dt) 
    end
    return k_capital
end

# Boundary error option 1
function boundary_condition_ergodic_1_step(p::Params, k_capital, C)
    @unpack A, α, δ, dt, k_min, k_max = p

    Y = A * k_capital .^ α # output
    I = Y - exp.(C(k_capital)) # investment
    dK_dt = I - δ .* k_capital # capital drift (also known as dK/dt)
    k_capital_t_plus_one = k_capital + (dt .* dK_dt)
    
    # we require that k_min < k(t+1) < k_max!
    error_lowerbound = max.(k_min .- k_capital_t_plus_one, 0)
    error_upperbound = max.(k_capital_t_plus_one .- k_max, 0)
    error = error_lowerbound .+ error_upperbound
    return error
end

# Boundary error option 2
function boundary_condition(p::Params, k_capital, C)
    @unpack A, α, δ, ρ, n_burned, k_min, epsilon = p
    
    Y = A * k_capital .^ α # output
    I = Y - exp.(C(k_capital)) # investment
    dK_dt = I - δ * k_capital # capital drift (also known as dK/dt)
    error = ifelse.(k_capital .< epsilon .&& dK_dt .< 0, dK_dt, 0)   

    return error
end

# Consumption error 
function c_error(p::Params, k_capital, V, C)
    @unpack γ = p

    v_prime = Zygote.gradient(k_capital -> sum(V(k_capital)), k_capital)[1] # dV/dk
    v_prime_max = max.(v_prime, 1e-7) 
    c_err = v_prime_max .^ (-1 ./ γ) - exp.(C(k_capital)) 

    return c_err
end

### LOSS FUNCTION

# Set global seed
Random.seed!(42)

function loss_function(V, C)

    # Training data
    k_capital = Flux.flatten(rand(Uniform(0.1, 10.0), 1000)) 

    ergodic_k_capital = ergodic_distribution(p, k_capital, C) 
    
    #error_v = hjb(p, ergodic_k_capital, V, C)
    #error_c = c_error(p, ergodic_k_capital, V, C)
    #error_b = boundary_condition(p, ergodic_k_capital, V, C)
    #error_b = boundary_condition_ergodic_1_step(p, ergodic_k_capital, C)

    #loss_v = mean(error_v)^2 # `mean` should do the same as `tf.reduce_mean` and `np.mean`
    #loss_c = mean(error_c)^2
    #loss_b = mean(error_b)^2
    #total_loss = loss_v + loss_c + loss_b

    # The following should result in fewer allocations(?)
    loss_v = mean(hjb(p, ergodic_k_capital, V, C))^2
    loss_c = mean(c_error(p, ergodic_k_capital, V, C))^2
    loss_b = mean(boundary_condition(p, ergodic_k_capital, C))^2
    # loss_b = mean(boundary_condition_ergodic_1_step(p, ergodic_k_capital, C))^2
    total_loss = loss_v + loss_c + loss_b

    return total_loss 
end

#### TRAIN THE MODEL

# Set up the neural network
model_vf = neural_net_vf(p)
model_cf = neural_net_cf(p)  

# Choose the appropriate optimizer
opt = ADAM() # can specify η learning rate and decay of momentum β::Tuple here

# Custom training loop
function custom_train!(loss_function, model_vf, model_cf, opt)

    # Define the trainable network parameters
    ps_vf = Flux.params(model_vf) # trainable parameters of model (VF)
    ps_cf = Flux.params(model_cf) # trainable parameters of model (CF)

    # This is for the total loss, we calculate the loss and gradient values using the pullback function
    total_loss_vf, grads_vf = Zygote.pullback(() -> loss_function(model_vf, model_cf), ps_vf)
    total_loss_cf, grads_cf = Zygote.pullback(() -> loss_function(model_vf, model_cf), ps_cf)

    # gradients in correct format for evaluation
    grads_v = grads_vf(one(total_loss_vf))      
    grads_c = grads_cf(one(total_loss_cf))

    # updates the parameters (same as optimizer.apply_gradients in TF) via ADAM
    Flux.Optimise.update!(opt, ps_vf, grads_v) 
    Flux.Optimise.update!(opt, ps_cf, grads_c)  
end

#### RUN THE MODEL

# Function to run model and track losses
function run_model(num_epochs)

    losses = []

    for i in 1:num_epochs
        custom_train!(loss_function, model_vf, model_cf, opt)
        push!(losses, loss_function(model_vf, model_cf))

        println("Epoch = $i")
    end
    return losses
end

# Run the model for certain amount of epochs
losses = run_model(10000) 

# Print the results
println("Sum of errors: ", losses[end])

#### PLOTTING
plot(losses, xaxis = :log, yaxis = :log) # Log-log scale

# Capital range
K = LinRange(0.1, 10.0, 1000);

# Value function plot
plot(K, model_vf(K')')

# Policy function plot
plot(K, model_cf(K')')