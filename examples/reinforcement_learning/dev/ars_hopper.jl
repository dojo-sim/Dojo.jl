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
vis=visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))

env = get_environment("hopper", vis=vis, dt = 0.05, damper = 2.0);
obs = reset(env)
render(env)

hp = HyperParameters(main_loop_size = 30, horizon = 0, n_directions = 6, b = 6, step_size = 0.02)

input_size = length(obs)
output_size = length(env.input_previous)
policy = Policy(input_size, output_size, hp)
normalizer = Normalizer(input_size)

train(env, policy, normalizer, hp)

traj = display_policy(env, policy, normalizer, hp)
traj = display_random_policy(env, hp)
visualize(env, traj)


# θbest005 = deepcopy(policy.θ)
# θbest = deepcopy(policy.θ)
# θgood = deepcopy(policy.θ)
# θbad = deepcopy(policy.θ)

# 5x = 1.1*1.5*2
# 2x = 1.4
