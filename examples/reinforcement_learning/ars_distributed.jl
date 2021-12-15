################################################################################
# Augmented Random Search ARS
################################################################################
using LinearAlgebra
using Statistics

import LinearAlgebra.normalize
import GeometryBasics.update

# ARS options: hyper parameters
@with_kw struct HyperParameters{T}
    main_loop_size::Int = 100
    horizon::Int = 200
    step_size::T = 0.02
    n_directions::Int = 16
    b::Int = 16
    noise::T = 0.03
    seed::Int = 1
    env_name::String = "halfcheetah"
end

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

function observe(normalizer::Normalizer{T}, x::AbstractVector{T}) where {T}
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

function Policy(input_size::Int, output_size::Int, hp::HyperParameters{T}; scale=1.0e-1) where {T}
    return Policy{T}(hp, scale * randn(output_size, input_size))
end

function evaluate(policy::Policy{T}, input::AbstractVector{T}) where {T}
    return policy.θ * input
end

function positive_perturbation(policy::Policy{T}, input::AbstractVector{T}, δ::AbstractMatrix{T}) where {T}
    return (policy.θ + policy.hp.noise .* δ) * input
end

function negative_perturbation(policy::Policy{T}, input::AbstractVector{T}, δ::AbstractMatrix{T}) where {T}
    return (policy.θ - policy.hp.noise .* δ) * input
end

function sample_δs(policy::Policy{T}) where {T}
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

function train(env::Environment, policy::Policy{T}, normalizer::Normalizer{T}, hp::HyperParameters{T}) where T
    envs = [deepcopy(env) for i = 1:Threads.nthreads()]
    for episode = 1:hp.main_loop_size
        # init deltas and rewards
        δs = sample_δs(policy)
        reward_positive = zeros(hp.n_directions)
        reward_negative = zeros(hp.n_directions)

        # positive directions
        Threads.@threads for k = 1:hp.n_directions
            seed(env, s = episode*k) #TODO I added this not sure if good or not
            state = reset(envs[Threads.threadid()])
            done = false
            num_plays = 0.
            while !done && num_plays < hp.horizon
                observe(normalizer, state)
                state = normalize(normalizer, state)
                action = positive_perturbation(policy, state, δs[k])
                state, reward, done, _ = step(envs[Threads.threadid()], action)
                reward = max(min(reward, 1), -1)
                reward_positive[k] += reward
                num_plays += 1
            end
        end

        # negative directions
        Threads.@threads for k = 1:hp.n_directions
            seed(env, s = episode*k) #TODO I added this not sure if good or not
            state = reset(envs[Threads.threadid()])
            done = false
            num_plays = 0.
            while !done && num_plays < hp.horizon
                observe(normalizer, state)
                state = normalize(normalizer, state)
                action = negative_perturbation(policy, state, δs[k])
                state, reward, done, _ = step(envs[Threads.threadid()], action)
                reward = max(min(reward, 1), -1)
                reward_negative[k] += reward
                num_plays += 1
            end
        end

        all_rewards = [reward_negative; reward_positive]
        σ_r = std(all_rewards)

        r_max = [max(reward_negative[k], reward_positive[k]) for k = 1:hp.n_directions]
        order = sortperm(r_max, rev = true)[1:hp.b]
        rollouts = [(reward_positive[k], reward_negative[k], δs[k]) for k = order]
        update(policy, rollouts, σ_r)
        
        reward_evaluation = mean(all_rewards)

        # finish, print:
        println("episode $episode reward_evaluation $reward_evaluation")
    end

    return nothing
end

