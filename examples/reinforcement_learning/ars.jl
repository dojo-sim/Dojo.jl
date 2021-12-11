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

function Policy(input_size::Int, output_size::Int, hp::HyperParameters{T}) where {T}
    return Policy{T}(hp, zeros(output_size, input_size))
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
    for episode = 1:hp.main_loop_size
        # init deltas and rewards
        δs = sample_δs(policy)
        reward_positive = zeros(hp.n_directions)
        reward_negative = zeros(hp.n_directions)

        # positive directions
        for k = 1:hp.n_directions
            state = reset(env)
            done = false
            num_plays = 0.
            while !done && num_plays < hp.horizon
                observe(normalizer, state)
                state = normalize(normalizer, state)
                action = positive_perturbation(policy, state, δs[k])
                state, reward, done, _ = step(env, action)
                reward = max(min(reward, 1), -1)
                reward_positive[k] += reward
                num_plays += 1
            end
        end

        # negative directions
        for k = 1:hp.n_directions
            state = reset(env)
            done = false
            num_plays = 0.
            while !done && num_plays < hp.horizon
                observe(normalizer, state)
                state = normalize(normalizer, state)
                action = negative_perturbation(policy, state, δs[k])
                state, reward, done, _ = step(env, action)
                reward = max(min(reward, 1), -1)
                reward_negative[k] += reward
                num_plays += 1
            end
        end

        all_rewards = [reward_negative; reward_positive]
        σ_r = std(all_rewards)

        # sort rollouts wrt max(r_pos, r_neg) and take (hp.b) best
        # scores = {k:max(r_pos, r_neg) for k,(r_pos,r_neg) in enumerate(zip(reward_positive,reward_negative))}
        # order = sorted(scores.keys(), key=lambda x:scores[x])[-hp.b:]
        # rollouts = [(reward_positive[k], reward_negative[k], deltas[k]) for k in order[::-1]]
        r_max = [max(reward_negative[k], reward_positive[k]) for k = 1:hp.n_directions]
        order = sortperm(r_max, rev = true)[1:hp.b]
        rollouts = [(reward_positive[k], reward_negative[k], δs[k]) for k = order]
        update(policy, rollouts, σ_r)
        @show scn.(policy.θ)

        # evaluate
        state = reset(env)
        done = false
        num_plays = 1.
        reward_evaluation = 0
        while !done && num_plays<hp.horizon
            observe(normalizer, state)
            state = normalize(normalizer, state)
            action = evaluate(policy, state)
            state, reward, done, _ = step(env, action)
            reward_evaluation += reward
            num_plays += 1
        end

        # finish, print:
        println("episode $episode reward_evaluation $reward_evaluation")
    end

    return nothing
end


# display learned policy
function display_policy(env::Environment, policy::Policy, normalizer::Normalizer, hp::HyperParameters)
    state = reset(env)
    done = false
    num_plays = 1.
    reward_evaluation = 0
    while !done && num_plays < hp.horizon
        render(env)
        sleep(env.mechanism.Δt)
        observe(normalizer, state)
        state = normalize(normalizer, state)
        action = evaluate(policy, state)
        state, reward, done, _ = step(env, action)
        reward_evaluation += reward
        num_plays += 1
    end
    close(env)
    return nothing
end
