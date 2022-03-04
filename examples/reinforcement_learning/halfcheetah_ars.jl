using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# ## Setup
using Dojo
using Random
using LinearAlgebra 
using JLD2
include(joinpath(@__DIR__, "algorithms/ars.jl")) # augmented random search

# ## Environment
env = get_environment(:halfcheetah, 
    timestep=0.05)
obs = reset(env)

# ## Open visualizer
open(env.vis)

# ## Training
train_times = Float64[]
rewards = Float64[]
policies = Matrix{Float64}[]
N = 5
for i = 1:N
    ## Reset environment
    env = get_environment(:halfcheetah, 
        timestep=0.05)
    obs = reset(env)

    ## Random policy
    hp = HyperParameters(
            main_loop_size=30, 
            horizon=80, 
            n_directions=6, 
            b=6, 
            step_size=0.02)
    input_size = length(obs)
    output_size = length(env.input_previous)
    normalizer = Normalizer(input_size)
    policy = Policy(input_size, output_size, hp)

    ## Train policy
    train_time = @elapsed train(env, policy, normalizer, hp)

    ## Evaluate policy
    reward = rollout_policy(policy.θ, env, normalizer, hp)

    # Cache
    push!(train_times, train_time)
    push!(rewards, reward)
    push!(policies, policy.θ)
end

## @save joinpath(@__DIR__, "results/halfcheetah_rl.jld2") train_times rewards policies
@load joinpath(@__DIR__, "results/halfcheetah_rl.jld2") train_times rewards policies

# Training statistics
N_best = 3
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

# ## Recover policy
hp = HyperParameters(
        main_loop_size=30, 
        horizon=80, 
        n_directions=6, 
        b=6, 
        step_size=0.02)
input_size = length(obs)
output_size = length(env.input_previous)
normalizer = Normalizer(input_size)

θ = policies_best[2]
traj = display_policy(env,
    ## policy,
    Policy(hp, θ),
    normalizer, hp)

for t = 1:length(traj)
    traj[t][2] += 2.0
end

# ## Visualizer policy
open(env.vis)
visualize(env, traj)
set_camera!(env.vis, 
    cam_pos=[0.0, -52.0, 0.0], 
    zoom=15)




