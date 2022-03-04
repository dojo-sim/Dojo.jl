using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# ## Setup
using Dojo
using Random
using LinearAlgebra 
using JLD2
include(joinpath(@__DIR__, "algorithms/ars.jl")) # augmented random search

# ## Ant
env = get_environment(:ant, 
    representation=:minimal, 
    gravity=-9.81, 
    timestep=0.05, 
    damper=50.0, 
    spring=25.0, 
    friction_coefficient=0.5,
    contact_feet=true, 
    contact_body=true)

obs = reset(env)
initialize!(env.mechanism, :ant,
    body_position=[0.0, 0.0, 1.0], 
    body_orientation=[0.0, 0.0, 0.0])
env.state .= get_minimal_state(env.mechanism)
render(env)

# ## Open visualizer
initialize!(env.mechanism, :ant)
open(env.vis)

# ## Set up policy
hp = HyperParameters(
        main_loop_size=100, 
        horizon=150, 
        n_directions=6, 
        b=6, 
        step_size=0.02)
input_size = length(obs)
output_size = length(env.input_previous)
normalizer = Normalizer(input_size)

# ## Training
train_times = Float64[]
rewards = Float64[]
policies = Matrix{Float64}[]
N = 5
for i = 1:N
    ## Reset environment
    env = get_environment(:ant, 
        representation=:minimal, 
        gravity=-9.81, 
        timestep=0.05, 
        damper=50.0, 
        spring=25.0, 
        friction_coefficient=0.5,
        contact_feet=true, 
        contact_body=true)
    obs = reset(env)

    ## Random policy
    Random.seed!(i)
    hp = HyperParameters(
        main_loop_size=100, 
        horizon=150, 
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

    ## Cache
    push!(train_times, train_time)
    push!(rewards, reward)
    push!(policies, policy.θ)
end

## @save joinpath(@__DIR__, "results/ant_rl.jld2") train_times rewards policies
@load joinpath(@__DIR__, "results/ant_rl.jld2") train_times rewards policies

# ## Training statistics
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

# ## Save/Load policy
## θ = policy.θ
## @save joinpath(@__DIR__, "ant_policy.jld2") θ
## @load joinpath(@__DIR__, "ant_policy.jld2") θ

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

# ## Visualize policy
## traj = display_random_policy(env, hp)
traj = display_policy(env, Policy(hp, θ), normalizer, hp)
visualize(env, traj)
open(env.vis)
