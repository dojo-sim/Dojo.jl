################################################################################
# Augmented Random Search ARS
################################################################################
using LinearAlgebra
using Statistics
using Folds

import LinearAlgebra.normalize
import Dojo.GeometryBasics.update

# ARS options: hyper parameters
Base.@kwdef struct HyperParameters{T}
    main_loop_size::Int = 100
    horizon::Int = 200
    step_size::T = 0.02
    n_directions::Int = 16
    b::Int = 16
    noise::T = 0.03
    seed::Int = 1
    env_name::String = "halfcheetah"
end
# assert self.b<=self.n_directions, "b must be <= n_directions"

# observation filter
mutable struct Normalizer{T}
    n::Vector{T}
    mean::Vector{T}
    mean_diff::Vector{T}
    var::Vector{T}
end

function Normalizer(num_inputs::Int)
    n = zeros(num_inputs)
    mean = zeros(num_inputs)
    mean_diff = zeros(num_inputs)
    var = zeros(num_inputs)
    return Normalizer{eltype(n)}(n, mean, mean_diff, var)
end

function observe(normalizer::Normalizer{T}, x::AbstractVector{T}) where T
    normalizer.n .+= 1
    last_mean = deepcopy(normalizer.mean)
    normalizer.mean .+= (x .- normalizer.mean) ./ normalizer.n
    normalizer.mean_diff .+= (x .- last_mean) .* (x .- normalizer.mean)
    normalizer.var .= max.(1e-2, normalizer.mean_diff ./ normalizer.n)
end

function normalize(normalizer, inputs)
    obs_mean = normalizer.mean
    obs_std = sqrt.(normalizer.var)
    return (inputs .- obs_mean) ./ obs_std
end

# linear policy
mutable struct Policy{T}
    hp::HyperParameters{T}
    θ::Matrix{T}
end

function Policy(input_size::Int, output_size::Int, hp::HyperParameters{T}; scale=1.0e-1) where T
    return Policy{T}(hp, scale * randn(output_size, input_size))
end

function evaluate(policy::Policy{T}, input::AbstractVector{T}) where T
    return policy.θ * input
end

function positive_perturbation(policy::Policy{T}, input::AbstractVector{T}, δ::AbstractMatrix{T}) where T
    return (policy.θ + policy.hp.noise .* δ) * input
end

function negative_perturbation(policy::Policy{T}, input::AbstractVector{T}, δ::AbstractMatrix{T}) where T
    return (policy.θ - policy.hp.noise .* δ) * input
end

function sample_δs(policy::Policy{T}) where T
    return [randn(size(policy.θ)) for i = 1:policy.hp.n_directions]
end

function update(policy::Policy, rollouts, σ_r)
    stp = zeros(size(policy.θ))
    for (r_pos, r_neg, d) in rollouts
        stp += (r_pos - r_neg) * d
    end
    policy.θ += policy.hp.step_size * stp ./ (σ_r * policy.hp.b)
    return nothing
end

function sample_policy(policy::Policy{T}) where T
    δ = [randn(size(policy.θ)) for i = 1:policy.hp.n_directions]
    θp = [policy.θ + policy.hp.noise .* δ[i] for i = 1:policy.hp.n_directions]
    θn = [policy.θ - policy.hp.noise .* δ[i] for i = 1:policy.hp.n_directions]
    return [θp..., θn...], δ
end

function eval_sample_policy(θ::Matrix{T}, input::AbstractVector{T}) where T
    return θ * input
end

function rollout_policy(θ::Matrix, env::Environment, normalizer::Normalizer, hp::HyperParameters;
        reset=reset)
    state = reset(env)
    rewards = 0.0
    done = false
    num_plays = 0
    while !done && num_plays < hp.horizon
        observe(normalizer, state)
        state = normalize(normalizer, state)
        action = eval_sample_policy(θ, state)
        state, reward, done, _ = step(env, action)
        reward = max(min(reward, 1), -1)
        rewards += reward
        num_plays += 1
    end
    return rewards
end

function train(env::Environment, policy::Policy{T}, normalizer::Normalizer{T},
        hp::HyperParameters{T}; distributed=false, usefolds=false, foldsexec=Folds.ThreadedEx(;basesize=1)) where T
    println("Training linear policy with Augmented Random Search (ARS)\n ")
    if distributed
        envs = [deepcopy(env) for i = 1:(2 * hp.n_directions)]
        normalizers = [deepcopy(normalizer) for i = 1:(2 * hp.n_directions)]
        hps = [deepcopy(hp) for i = 1:(2 * hp.n_directions)]
        print("  $(nprocs()) processors")
    elseif usefolds
        envs = [deepcopy(env) for i = 1:(2*hp.n_directions)]
        print("  $(Threads.nthreads()) threads with Folds")
    else
        envs = [deepcopy(env) for i = 1:(Threads.nthreads())]
        print(" $(Threads.nthreads()) ")
    end

    # pre-allocate for rewards
    rewards = zeros(2 * hp.n_directions)

    for episode = 1:hp.main_loop_size
        # initialize deltas and rewards
        θs, δs = sample_policy(policy)

        # evaluate policies
        roll_time = @elapsed begin
            if distributed
                rewards .= pmap(rollout_policy, θs, envs, normalizers, hps)
            elseif usefolds
                @assert length(envs) == size(θs, 1) "$(length(envs))"
                Folds.map!(rewards, θs, envs, foldsexec) do θ, env
                    rollout_policy(θ, env, normalizer, hp)
                end
            else
                Threads.@threads for k = 1:(2 * hp.n_directions)
                    rewards[k] = rollout_policy(θs[k], envs[Threads.threadid()], normalizer, hp)
                end
            end
        end
        # reward evaluation
        r_max = [max(rewards[k], rewards[hp.n_directions + k]) for k = 1:hp.n_directions]
        σ_r = std(rewards)
        order = sortperm(r_max, rev = true)[1:hp.b]
        rollouts = [(rewards[k], rewards[hp.n_directions + k], δs[k]) for k = order]

        # policy update
        update(policy, rollouts, σ_r)

        # finish, print:
        println("episode $episode reward_evaluation $(mean(rewards)). Took $(roll_time) seconds")
    end

    return nothing
end

# display learned policy
function display_policy(env::Environment, policy::Policy, normalizer::Normalizer, hp::HyperParameters; rendering = false, reset=reset)
    obs = reset(env)
    traj = [copy(env.state)]
    done = false
    num_plays = 1.
    reward_evaluation = 0
    while !done && num_plays < hp.horizon
        rendering && render(env)
        sleep(env.mechanism.timestep)
        observe(normalizer, obs)
        push!(traj, copy(env.state))
        obs = normalize(normalizer, obs)
        action = evaluate(policy, obs)
        obs, reward, done, _ = step(env, action)
        reward_evaluation += reward
        num_plays += 1
    end
    close(env)

    return traj
end

# display learned policy
function display_random_policy(env::Environment, hp::HyperParameters; rendering = false)
    obs = reset(env)
    traj = [copy(env.state)]
    done = false
    num_plays = 1.
    reward_evaluation = 0
    while !done && num_plays < hp.horizon
        rendering && render(env)
        sleep(env.mechanism.timestep)
        push!(traj, copy(env.state))
        action = sample(env.input_space)
        obs, reward, done, _ = step(env, action)
        reward_evaluation += reward
        num_plays += 1
    end
    close(env)

    return traj
end
