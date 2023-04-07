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
    # Rollout data
    X = [zeros(nx) for k = 1:hp.horizon+1]
    A = [zeros(nu) for k = 1:hp.horizon]
    fx = [zeros(nx,nx) for k = 1:hp.horizon]
    fu = [zeros(nx,nu) for k = 1:hp.horizon]
    # Gradient data
    ∇θ = zeros(1,nθ)
    ∂x∂θ = [zeros(nx, nθ) for k = 1:hp.horizon]

    for episode = 1:hp.main_loop_size
        # Stem
        state = reset(env)
        X[1] .= copy(state)
        done = false
        sum_reward = 0
        num_plays = 0
        for k = 1:hp.horizon
            done && break
            observe(normalizer, state)
            state = normalize(normalizer, state)
            action = evaluate(policy, state)
            A[k] .= copy(action)
            state, reward, done, _ = step(env, action, diff = true)
            X[k+1] .= copy(state)
            fx[k] .= copy(env.dynamics_jacobian_state)
            fu[k] .= copy(env.dynamics_jacobian_input)
            reward = max(min(reward, 1), -1)
            sum_reward += reward
            num_plays += 1
        end

        ∇θ .= 0.0
        ∂u∂θ = [cat([X[k]' for i = 1:output_size]..., dims = (1,2)) for k = 1:num_plays]
        ∂r∂x = [FiniteDiff.finite_difference_jacobian(x -> [-cost(env, x, A[k])], X[k]) for k = 1:num_plays]
        ∂r∂u = [FiniteDiff.finite_difference_jacobian(u -> [-cost(env, X[k], u)], A[k]) for k = 1:num_plays]
        for k = 2:num_plays
            ∂x∂θ[k] .= fx[k-1] * ∂x∂θ[k-1] + fu[k-1] * ∂u∂θ[k-1] #TODO
        end
        for k = 1:num_plays
            if norm(∇θ, Inf) > 8e2
                println("grad up to step $k")
                break
            end
            ∇θ += ∂r∂x[k] * ∂x∂θ[k] + ∂r∂u[k] * ∂u∂θ[k]
        end
        ∇ = transpose(reshape(∇θ, (input_size, output_size)))
        (norm(∇, Inf) < 1e3) && gradient_update(policy, ∇)

        # finish, print:
        println("episode $episode ∇∞ $(scn(norm(∇, Inf))) r $sum_reward")
    end
    return nothing
end

function gradient_update(policy::Policy, ∇)
    policy.θ += policy.hp.step_size * ∇ ./ norm(∇, Inf)
    return nothing
end
