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
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))

# env = make("ant", vis = vis, mode=:min, g=-9.81, dt=0.05, damper=5.0, spring=5.0, cf = 0.5,
env = make("ant", vis = vis, mode=:min, g=-9.81, dt=0.05, damper=50.0, spring=30.0, cf = 0.5,
    contact=true, contact_body=true)
obs = reset(env)
initializeant!(env.mechanism, pos = [1.3,0,0], rot = [0,0,0.])
env.x .= getMinState(env.mechanism)
render(env)

hp = HyperParameters(main_loop_size = 30, horizon = 150, n_directions = 6, b = 6, step_size = 0.02)

input_size = length(obs)
output_size = length(env.u_prev)
# policy = Policy(input_size, output_size, hp)
# normalizer = Normalizer(input_size)
train(env, policy, normalizer, hp)



traj = display_policy(env, policy, normalizer, hp)
visualize(env, traj)
# traj = display_random_policy(env, hp)
# visualize(env, traj)
