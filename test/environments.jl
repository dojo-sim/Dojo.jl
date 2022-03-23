environments = [
    :ant,
    :atlas,
    :block,
    :block2D,
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
    :block2D,
    :hopper,
    :panda,
    :rexhopper,
    :walker,
]

@testset "$name" for name in environments
    env = get_environment(name)
    if !(name in throw_envs)
        @test size(reset(env)) == (env.observation_space.n,)
        o, r, d, i = step(env, Dojo.sample(env.input_space))
        @test size(o) == (env.observation_space.n,)
    end
end

using Test
for name in environments
    @show name
    env = get_environment(name)
    if !(name in throw_envs)
        @test size(reset(env)) == (env.observation_space.n,)
        o, r, d, i = step(env, Dojo.sample(env.input_space))
        @test size(o) == (env.observation_space.n,)
    end
end

env = get_environment(:ant)
env = get_environment(:atlas)
env = get_environment(:block)
env = get_environment(:block2D)
env = get_environment(:cartpole)
env = get_environment(:halfcheetah)
env = get_environment(:hopper)
env = get_environment(:panda)
env = get_environment(:pendulum)
env = get_environment(:quadruped)
env = get_environment(:raiberthopper)
env = get_environment(:rexhopper)
env = get_environment(:walker)
@test size(reset(env)) == (env.observation_space.n,)
o, r, d, i = step(env, Dojo.sample(env.input_space))
@test size(o) == (env.observation_space.n,)

:atlas
:block
