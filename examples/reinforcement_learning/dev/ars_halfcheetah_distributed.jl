using Distributed 
@show nprocs()
# addprocs(3)

@everywhere begin
    # Utils
    function module_dir()
        return joinpath(@__DIR__, "..", "..", "..")
    end

    # Activate package
    using Pkg
    Pkg.activate(module_dir())
    using Dojo
end 

env = make("halfcheetah", dt=0.05)
obs = reset(env)
hp = HyperParameters(main_loop_size=30, horizon=80, n_directions=6, b=6, step_size=0.02)
input_size = length(obs)
output_size = length(env.u_prev)
policy = Policy(input_size, output_size, hp)
normalizer = Normalizer(input_size)

train(env, policy, normalizer, hp, distributed=true)

# Visualizer
vis = Visualizer()
open(env.vis)
traj = display_policy(env, policy, normalizer, hp)
visualize(env, traj)