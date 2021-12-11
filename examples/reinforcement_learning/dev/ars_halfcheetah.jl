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



env = make("pendulum", g = -1.00, damper = 0.4, max_torque = 12.0, vis = vis)
obs = reset(env)
render(env)
input_size = length(obs)
output_size = length(env.u_prev)
hp = Hp13(main_loop_size = 200, horizon = 100, n_directions = 4, b = 4, step_size = 0.05)
normalizer = Normalizer13(input_size)
policy = Policy13(input_size, output_size, hp)
train(env, policy, normalizer, hp)
display_policy(env, policy, normalizer, hp)

# policy.θ = 1.5*[1.0 -0.3]
