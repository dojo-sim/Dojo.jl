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

# ## Reset environment
env = get_environment(:halfcheetah, timestep=0.05)
obs = reset(env)

# ## Random policy
Random.seed!(38)
hp = HyperParameters(
    main_loop_size=10,
    horizon=80,
    n_directions=6,
    b=6,
    step_size=0.02)
input_size = length(obs)
output_size = length(env.input_previous)
policy = Policy(input_size, output_size, hp)
normalizer = Normalizer(input_size)
reward = rollout_policy(policy.θ, env, normalizer, hp, reset=env->reset(env,reset_noise_scale=0.0))
# ## Train policy
train_time = @elapsed train(env, policy, normalizer, hp)

# ## Evaluate policy
reward = rollout_policy(policy.θ, env, normalizer, hp)

# @save joinpath(@__DIR__, "results/ars_demonstration.jld2") train_time reward policy normalizer
# @load joinpath(@__DIR__, "results/ars_demonstration.jld2") train_time reward policy normalizer

# ## Visualizer policy
open(env.vis)
visualize(env, traj)
set_camera!(env.vis,
    cam_pos=[0.0, -52.0, 0.0],
    zoom=1)

convert_frames_to_video_and_gif("halfcheetah_square")
