# ## Setup
using Dojo
using DojoEnvironments
using Random
using LinearAlgebra 
using JLD2

using Statistics

import LinearAlgebra.normalize

# ## Parameters
rng = MersenneTwister(1)

Base.@kwdef struct HyperParameters{T}
    main_loop_size::Int = 100
    horizon::Int = 200
    step_size::T = 0.05
    n_directions::Int = 16
    b::Int = 16
    noise::T = 0.1
end

# ## Observation filter
mutable struct Normalizer{T}
    n::Vector{T}
    mean::Vector{T}
    mean_diff::Vector{T}
    var::Vector{T}

    function Normalizer(num_inputs::Int)
        n = zeros(num_inputs)
        mean = zeros(num_inputs)
        mean_diff = zeros(num_inputs)
        var = zeros(num_inputs)
        new{eltype(n)}(n, mean, mean_diff, var)
    end
end

function normalize(normalizer, inputs)
    obs_std = sqrt.(normalizer.var)
    return (inputs .- normalizer.mean) ./ obs_std
end

function observe!(normalizer, inputs)
    normalizer.n .+= 1
    last_mean = deepcopy(normalizer.mean)
    normalizer.mean .+= (inputs .- normalizer.mean) ./ normalizer.n
    normalizer.mean_diff .+= (inputs .- last_mean) .* (inputs .- normalizer.mean)
    normalizer.var .= max.(1e-2, normalizer.mean_diff ./ normalizer.n)
end


# ## Mechanism
env = get_environment(:ant_ars;
    horizon=100,
    gravity=-9.81, 
    timestep=0.05, 
    dampers=50.0, 
    springs=25.0, 
    friction_coefficient=0.5,
    contact_feet=true, 
    contact_body=true);

# ## Reset and rollout functions
function reset_state!(env)
    initialize!(env.mechanism, :ant)
    return
end

function rollout_policy(θ::Matrix, env, normalizer::Normalizer, parameters::HyperParameters; record=false)
    reset_state!(env)
    rewards = 0.0
    for k=1:parameters.horizon
        state = DojoEnvironments.get_state(env)
        x = state[1:28] # minimal state without contacts
        observe!(normalizer, state)
        state = normalize(normalizer, state)
        action = θ * state

        
        step!(env, x, action; record, k)
        state_after = DojoEnvironments.get_state(env)
        x_pos_before = x[1]
        x_pos_after = state_after[1]
        forward_reward = 100 * (x_pos_after - x_pos_before) / env.mechanism.timestep
        control_cost = (0.05/10 * action' * action)[1]
        contact_cost = 0.0
        for contact in env.mechanism.contacts
            contact_cost += 0.5 * 1.0e-3 * max(-1, min(1, contact.impulses[2][1]))^2.0
        end
        survive_reward = 0.05

        reward = forward_reward - control_cost - contact_cost + survive_reward
        # reward = max(min(reward, 1), -1)
        rewards += reward

        if !(all(isfinite.(state_after)) && (state_after[3] >= 0.2) && (state_after[3] <= 1))
            println("  failed")
            break
        end
    end

    println("  x_pos: "*string(DojoEnvironments.get_state(env)[1]))

    return rewards
end


# ## Training functions
mutable struct Policy{T}
    hp::HyperParameters{T}
    θ::Matrix{T}

    function Policy(input_size::Int, output_size::Int, hp::HyperParameters{T}; scale=0.2) where T
        new{T}(hp, scale * randn(output_size, input_size))
    end
end

function sample_policy(policy::Policy{T}) where T
    δ = [randn(size(policy.θ)) for i = 1:policy.hp.n_directions]
    θp = [policy.θ + policy.hp.noise .* δ[i] for i = 1:policy.hp.n_directions]
    θn = [policy.θ - policy.hp.noise .* δ[i] for i = 1:policy.hp.n_directions]
    return [θp..., θn...], δ
end

function update!(policy::Policy, rollouts, σ_r)
    stp = zeros(size(policy.θ))
    for (r_pos, r_neg, d) in rollouts
        stp += (r_pos - r_neg) * d
    end
    policy.θ += policy.hp.step_size * stp ./ (σ_r * policy.hp.b)
    return
end

function train(env, policy::Policy{T}, normalizer::Normalizer{T}, hp::HyperParameters{T}) where T

    println("Training linear policy with Augmented Random Search (ARS)\n ")

    # pre-allocate for rewards
    rewards = zeros(2 * hp.n_directions)

    for episode = 1:hp.main_loop_size
        # initialize deltas and rewards
        θs, δs = sample_policy(policy)

        # evaluate policies
        roll_time = @elapsed begin
            for k = 1:(2 * hp.n_directions)
                rewards[k] = rollout_policy(θs[k], env, normalizer, hp)
            end
        end
        # reward evaluation
        r_max = [max(rewards[k], rewards[hp.n_directions + k]) for k = 1:hp.n_directions]
        σ_r = std(rewards)
        order = sortperm(r_max, rev = true)[1:hp.b]
        rollouts = [(rewards[k], rewards[hp.n_directions + k], δs[k]) for k = order]

        # policy update
        update!(policy, rollouts, σ_r)

        # finish, print:
        println("episode $episode reward_evaluation $(mean(rewards)). Took $(roll_time) seconds")
        rollout_policy(policy.θ, env, normalizer, hp; record=true)
        visualize(env)
    end

    return nothing
end

# ## Training
train_times = Float64[]
rewards = Float64[]
policies = Matrix{Float64}[]
N = 1
for i = 1:N
    ## Random policy
    hp = HyperParameters(
        main_loop_size=20, 
        horizon=100, 
        n_directions=6, 
        b=6, 
        step_size=0.05)
    input_size = 37
    output_size = 8
    normalizer = Normalizer(input_size)
    policy = Policy(input_size, output_size, hp)

    ## Train policy
    train_time = @elapsed train(env, policy, normalizer, hp)

    ## Evaluate policy
    reward = rollout_policy(policy.θ, env, normalizer, hp)

    ## Cache
    push!(train_times, train_time)
    push!(rewards, reward)
    push!(policies, policy.θ)
end


# ## Training statistics
N_best = 1
max_idx = sortperm(rewards, 
    lt=Base.isgreater)
train_time_best = (train_times[max_idx])[1:N_best]
rewards_best = (rewards[max_idx])[1:N_best]
policies_best = (policies[max_idx])[1:N_best]

@show rewards
@show mean(train_time_best)
@show std(train_time_best)
@show mean(rewards)
@show std(rewards)

# ## Save/Load policy
# θ = policy.θ
# @save joinpath(@__DIR__, "ant_policy.jld2") θ
# @load joinpath(@__DIR__, "ant_policy.jld2") θ

# ## Recover policy
hp = HyperParameters(
        main_loop_size=20, 
        horizon=100, 
        n_directions=6, 
        b=6, 
        step_size=0.05)
input_size = 37
output_size = 8
normalizer = Normalizer(input_size)
θ = policies_best[1] 


function controller!(mechanism, k)
    state = DojoEnvironments.get_state(env)
    observe!(normalizer, state)
    state = normalize(normalizer, state)
    action = θ * state
    set_input!(mechanism, [zeros(6);action])
end

# # ## Visualize policy
reset_state!(env)
storage = simulate!(env.mechanism, 5, controller!; record=true)
rollout_policy(θ, env, normalizer, hp; record=true)
visualize(env.mechanism, storage)