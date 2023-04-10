environments = [
    :ant_ars,
    :cartpole_dqn,
    :quadruped_sampling,
]

for name in environments 
    env = get_environment(name) 
    x = get_state(env)
    u = input_map(env, nothing)
    step!(env, x)
    step!(env, x, u)
    simulate!(env, 0.5; record=true)
    visualize(env)
    @test true
end
