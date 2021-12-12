# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))


include("../ags.jl")







function train(env::Environment, policy::Policy{T}, normalizer::Normalizer{T}, hp::HyperParameters{T}) where T
    envs = [deepcopy(env) for i = 1:Threads.nthreads()]

    output_size, input_size = size(policy.θ)
    nx = input_size
    nu = output_size
    nθ = output_size * input_size


    for episode = 1:hp.main_loop_size
        # init deltas and rewards
        # δs = sample_δs(policy)
        # reward_positive = zeros(hp.n_directions)
        # reward_negative = zeros(hp.n_directions)

        # Stem
        seed(env, s = episode)
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
            ∂x∂θ[k] = fx[k-1] * ∂x∂θ[k-1] + fu[k-1] * ∂u∂θ[k-1]
        end
        for k = 1:num_plays
            ∇θ += ∂r∂x[k] * ∂x∂θ[k] + ∂r∂u[k] * ∂u∂θ[k]
        end
        ∇ = transpose(reshape(∇θ, (input_size, output_size)))
        @show mean(abs.(∇))
        gradient_update(policy, ∇)

        # all_rewards = [reward_negative; reward_positive]
        # σ_r = std(all_rewards)
        #
        # # sort rollouts wrt max(r_pos, r_neg) and take (hp.b) best
        # # scores = {k:max(r_pos, r_neg) for k,(r_pos,r_neg) in enumerate(zip(reward_positive,reward_negative))}
        # # order = sorted(scores.keys(), key=lambda x:scores[x])[-hp.b:]
        # # rollouts = [(reward_positive[k], reward_negative[k], deltas[k]) for k in order[::-1]]
        # r_max = [max(reward_negative[k], reward_positive[k]) for k = 1:hp.n_directions]
        # order = sortperm(r_max, rev = true)[1:hp.b]
        # rollouts = [(reward_positive[k], reward_negative[k], δs[k]) for k = order]
        # update(policy, rollouts, σ_r)
        # # @show scn.(policy.θ)

        # finish, print:
        println("episode $episode ∇∞ $(scn(norm(∇, Inf))) r $reward_evaluation")
    end
    return nothing
end

function gradient_update(policy::Policy, ∇)
    # @show norm(∇, 1) / (18*6)
    # @show mean(abs.(∇))
    # @show norm(∇, Inf)
    # policy.θ += policy.hp.step_size * ∇ #./ σ_r
    policy.θ += policy.hp.step_size * ∇ ./ norm(∇, Inf)
    return nothing
end

policy = Policy(input_size, output_size, hp)
normalizer = Normalizer(input_size)


env = make("halfcheetah", vis = vis, dt = 0.05)
obs = reset(env)
render(env)
input_size = length(obs)
output_size = length(env.u_prev)
hp = HyperParameters(main_loop_size = 30, horizon = 80, n_directions = 6, b = 6, step_size = 0.005)
train(env, policy, normalizer, hp)
traj = display_policy(env, policy, normalizer, hp)
visualize(env, traj)




# θbest005 = deepcopy(policy.θ)
# θbest = deepcopy(policy.θ)
# θgood = deepcopy(policy.θ)
# θbad = deepcopy(policy.θ)

# 5x = 1.1*1.5*2
# 2x = 1.4



# x0 = rand(18)
# u0 = rand(6)
# FiniteDiff.finite_difference_gradient(x -> cost(env, x, u), x0)
# FiniteDiff.finite_difference_gradient(u -> cost(env, x, u), u0)
#
# cat([rand(4)' for i = 1:2]..., dims = (1,2))
#
a = [1 2 3; 4 5 6]
ra = reshape(vec(a'), (1,6))
transpose(reshape(vec(a'), (3,2)))
