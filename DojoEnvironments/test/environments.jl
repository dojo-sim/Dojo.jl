environments = [
    :ant_ars,
    :cartpole_dqn,
    :quadruped_sampling,
]

for name in environments 
    env = get_environment(name; horizon=2) 
    x = get_state(env)
    u = DojoEnvironments.input_map(env, nothing)
    step!(env, x)
    step!(env, x, u)
    simulate!(env; record=true)
    visualize(env)
    @test true
end
