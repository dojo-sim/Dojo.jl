using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# ## Setup
using Dojo
using Random
using LinearAlgebra
using JLD2
include(joinpath(@__DIR__, "algorithms/ars.jl")) # augmented random search

# # Training

# ## Environment
env = get_environment(:halfcheetah, timestep=0.05)
obs = reset(env);

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
normalizer = Normalizer(input_size);

# ## Train policy
train_time = @elapsed train(env, policy, normalizer, hp)

# ## Evaluate policy
reward = rollout_policy(policy.θ, env, normalizer, hp, reset=env->reset(env,reset_noise_scale=0.0))

# @save joinpath(@__DIR__, "results/ars_demonstration_dynamic.jld2") train_time reward policy normalizer
@load joinpath(@__DIR__, "results/ars_demonstration_dynamic.jld2") train_time reward policy normalizer

# ## rollout trained policy
traj = display_policy(env,
    policy,
    normalizer,
    hp,
    reset=env->reset(env,reset_noise_scale=0.0))

for t = 1:length(traj)
    traj[t][2] += -0.50
end


# ## Open Visualizer
open(env.vis)

# ## Visualize Policy
visualize(env, traj)
set_camera!(env.vis,
    cam_pos=[0.0, -35.0, 0.0],
    zoom=15)

# set_camera!(env.vis,
#     cam_pos=[3.0, -6.0, 2.0],
#     zoom=2)
# render_static(env.vis)
# open(joinpath(@__DIR__, "results/ars_demonstration_dynamic.html"), "w") do file
#     write(file, static_html(env.vis))
# end
include("../generate_notebooks.jl")
