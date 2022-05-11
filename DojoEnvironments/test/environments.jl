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
        o, r, d, i = step(env, clamp.(DojoEnvironments.sample(env.input_space), -0.1, 0.1))
        @test size(o) == (env.observation_space.n,)
    end
end
