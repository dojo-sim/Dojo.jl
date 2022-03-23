environments = [
    :ant,
    :atlas,
    :block,
    :block2d,
    :cartpole,
    :halfcheetah,
    :hopper,
    :panda,
    :pendulum,
    :quadruped,
    :raiberthopper,
    :rexhopper,
    :walker,
]

throw_envs = [
    ]

@testset "$name" for name in environments
    env = get_environment(name)
    if !(name in throw_envs)
        @test size(reset(env)) == (env.observation_space.n,)
        o, r, d, i = step(env, clamp.(Dojo.sample(env.input_space), -0.1, 0.1))
        @test size(o) == (env.observation_space.n,)
    end
end
#
# using Test
# for name in environments
#     @show name
#     env = get_environment(name)
#     if !(name in throw_envs)
#         @test size(reset(env)) == (env.observation_space.n,)
#         o, r, d, i = step(env, clamp.(Dojo.sample(env.input_space), -0.1, 0.1))
#         @test size(o) == (env.observation_space.n,)
#     end
# end
#
# env.input_space
#
# env = get_environment(:ant)
# env = get_environment(:atlas)
# env = get_environment(:block)
# env = get_environment(:block2d)
# env = get_environment(:cartpole)
# env = get_environment(:halfcheetah)
# env = get_environment(:hopper)
# env = get_environment(:panda)
# env = get_environment(:pendulum)
# env = get_environment(:quadruped)
# env = get_environment(:raiberthopper)
# env = get_environment(:rexhopper)
# env = get_environment(:walker)
#
# reset(env)
# @test size(reset(env)) == (env.observation_space.n,)
# o, r, d, i = step(env, Dojo.sample(env.input_space))
# @test size(o) == (env.observation_space.n,)
#
# :atlas
# :block
# type2symbol(Block2D)
#
# vis = Visualizer()
# open(vis)
#
# initialize!(env.mechanism, :block)
# type2symbol(typeof(env))
# storage = simulate!(env.mechanism, 1.0)
# visualize(env.mechanism, storage, vis=vis)
#
# clamp.(rand(12), -0.01, 0.01)
#
# for name in [
#     :block2d,
#     :hopper,
#     :panda,
#     :rexhopper,
#     :walker,
# ]
#     @show name
#     env = get_environment(name)
#     @test size(reset(env)) == (env.observation_space.n,)
#     o, r, d, i = step(env, clamp.(Dojo.sample(env.input_space), -0.1, 0.1))
#     @test size(o) == (env.observation_space.n,)
# end
#
# env = get_environment(:rexhopper, infeasible_control=false)
# env = get_environment(:quadruped, infeasible_control=false)
# env.control_map
# input_dimension(env.mechanism)
# input_dimension(env.mechanism)
#
# env.num_inputs
#
#
# get_control_mask(14, [4,5,6,7,9])
#
# [zeros(nu, 3) I(4)]
