# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using MeshCat

# Open visualizer
vis=visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))

env = get_environment("halfcheetah", vis=vis, dt = 0.05)
obs = reset(env)
render(env)

hp = HyperParameters(main_loop_size = 30, horizon = 200, n_directions = 6, b = 6, step_size = 0.02)

input_size = length(obs)
output_size = length(env.input_previous)
policy = Policy(input_size, output_size, hp)
normalizer = Normalizer(input_size)

train(env, policy, normalizer, hp)
traj = display_policy(env, policy, normalizer, hp)

visualize(env, traj)
