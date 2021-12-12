################################################################################
# Augmented Gradient Search ARS
################################################################################
using LinearAlgebra
using Statistics

import LinearAlgebra.normalize
import GeometryBasics.update

function train(env::Environment, policy::Policy{T}, normalizer::Normalizer{T}, hp::HyperParameters{T}) where T
    envs = [deepcopy(env) for i = 1:Threads.nthreads()]

    output_size, input_size = size(policy.θ)
    nx = input_size
    nu = output_size
    nθ = output_size * input_size

    for episode = 1:hp.main_loop_size
        # Stem
        state = reset(env)
        X = [copy(state),]
        A = []
        fx = []
        fu = []
        done = false
        reward_evaluation = 0
        num_plays = 0
        while !done && num_plays < hp.horizon
            observe(normalizer, state)
            state = normalize(normalizer, state)
            action = evaluate(policy, state)
            push!(A, action)
            state, reward, done, _ = step(env, action, diff = true)
            push!(X, copy(state))
            push!(fx, copy(env.fx))
            push!(fu, copy(env.fu))
            reward = max(min(reward, 1), -1)
            reward_evaluation += reward
            num_plays += 1
        end

        ∇θ = zeros(1,nθ)
        ∂x∂θ = [zeros(nx, nθ) for k = 1:num_plays]
        ∂u∂θ = [cat([X[k]' for i = 1:output_size]..., dims = (1,2)) for k = 1:num_plays]
        ∂r∂x = [FiniteDiff.finite_difference_jacobian(x -> [-cost(env, x, A[k])], X[k]) for k = 1:num_plays]
        ∂r∂u = [FiniteDiff.finite_difference_jacobian(u -> [-cost(env, X[k], u)], A[k]) for k = 1:num_plays]
        for k = 2:num_plays
            ∂x∂θ[k] = fx[k-1] * ∂x∂θ[k-1] + fu[k-1] * ∂u∂θ[k-1] #TODO
        end
        for k = 1:num_plays
            ∇θ += ∂r∂x[k] * ∂x∂θ[k] + ∂r∂u[k] * ∂u∂θ[k]
        end
        ∇ = transpose(reshape(∇θ, (input_size, output_size)))
        (norm(∇, Inf) < 1e3) && gradient_update(policy, ∇)

        # finish, print:
        println("episode $episode ∇∞ $(scn(norm(∇, Inf))) r $reward_evaluation")
    end
    return nothing
end

function gradient_update(policy::Policy, ∇)
    @show "update"
    policy.θ += policy.hp.step_size * ∇ ./ norm(∇, Inf)
    return nothing
end
