# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))
# include("../ars_distributed.jl")

env = make("halfcheetah", dt=0.05)
obs = reset(env)

hp = HyperParameters(main_loop_size=1, horizon=10, n_directions=1, b=1, step_size=0.02)

input_size = length(obs)
output_size = length(env.u_prev)
policy = Policy(input_size, output_size, hp)
normalizer = Normalizer(input_size)

train(env, policy, normalizer, hp)

