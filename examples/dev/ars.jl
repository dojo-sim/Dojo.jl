################################################################################
# Augmented Random Search ARS
################################################################################

import LinearAlgebra.normalize
import GeometryBasics.update
using LinearAlgebra
using Statistics


# ARS options: hyper parameters
@with_kw struct Hp13{T}
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
mutable struct Normalizer13{T}
    n::Vector{T}
    mean::Vector{T}
    mean_diff::Vector{T}
    var::Vector{T}
end

function Normalizer13(num_inputs::Int)
    n = zeros(num_inputs)
    mean = zeros(num_inputs)
    mean_diff = zeros(num_inputs)
    var = zeros(num_inputs)
    return Normalizer13{eltype(n)}(n, mean, mean_diff, var)
end

function observe(normalizer::Normalizer13{T}, x::AbstractVector{T}) where {T}
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
mutable struct Policy13{T}
    hp::Hp13{T}
    θ::Matrix{T}
end

function Policy13(input_size::Int, output_size::Int, hp::Hp13{T}) where {T}
    return Policy13{T}(hp, zeros(output_size, input_size))
end

function evaluate(policy::Policy13{T}, input::AbstractVector{T}) where {T}
    return policy.θ * input
end

function positive_perturbation(policy::Policy13{T}, input::AbstractVector{T}, δ::AbstractMatrix{T}) where {T}
    return (policy.θ + policy.hp.noise .* δ) * input
end

function negative_perturbation(policy::Policy13{T}, input::AbstractVector{T}, δ::AbstractMatrix{T}) where {T}
    return (policy.θ - policy.hp.noise .* δ) * input
end

function sample_δs(policy::Policy13{T}) where {T}
    return [randn(size(policy.θ)) for i = 1:policy.hp.n_directions]
end

function update(policy::Policy13, rollouts, σ_r)
    stp = zeros(size(policy.θ))
    for (r_pos, r_neg, d) in rollouts
        stp += (r_pos - r_neg) * d
    end
    policy.θ += policy.hp.step_size * stp ./ (σ_r * policy.hp.b)
    return nothing
end

function train(env::Environment, policy::Policy13{T}, normalizer::Normalizer13{T}, hp::Hp13{T}) where T
    for episode = 1:hp.main_loop_size
        @show episode
        # init deltas and rewards
        δs = sample_δs(policy)
        reward_positive = zeros(hp.n_directions)
        reward_negative = zeros(hp.n_directions)

        # positive directions
        for k = 1:hp.n_directions
            # println("positive direction: k ", k)
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
            # @show reward_positive[k]
        end

        # negative directions
        for k = 1:hp.n_directions
            # println("negative direction: k ", k)
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
        # @show scn(rollouts[1][1]), scn(rollouts[1][2]), scn.(rollouts[1][3])
        @show scn.(policy.θ)
        # update policy:
        update(policy, rollouts, σ_r)

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
function display_policy(env::Environment, policy::Policy13, normalizer::Normalizer13, hp::Hp13)
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

env0 = make("pendulum", g = -1.00, damper = 0.4, max_torque = 20.0)
open(env0.vis)
obs0 = reset(env0)
render(env0)
input_size0 = length(obs0)
output_size0 = length(env0.u_prev)
hp0 = Hp13(main_loop_size = 200, horizon = 100, n_directions = 4, b = 4, step_size = 0.1)
normalizer0 = Normalizer13(input_size0)
policy0 = Policy13(input_size0, output_size0, hp0)
train(env0, policy0, normalizer0, hp0)
display_policy(env0, policy0, normalizer0, hp0)
policy0.θ * _get_obs(env)

# policy0.θ = 1.5*[1.0 -0.3]

obs0 = reset(env0, x = [-0.1,0.])
render(env0)
cost(env0, [-0.1, 0.], [0.])

policy0.θ.= []
cost(env, [0, 0.0], [0.0])


angle_normalize( 0.1 * π)
angle_normalize( 1.1 * π)
angle_normalize(-0.9 * π)
plot(-3:0.01:3, angle_normalize.(-pi .+ pi .* (-3:0.01:3)).^2)
