environments = [
    :ant_ars,
    :cartpole_dqn,
    :pendulum,
    :quadruped_waypoint
    :quadruped_sampling,
    :quadrotor_waypoint,
    :uuv_waypoint,
    :youbot_waypoint,
]

for name in environments 
    env = get_environment(name; horizon=2)
    x = get_state(env)
    u = DojoEnvironments.input_map(env, nothing)
    set_input!(env, nothing)
    step!(env, x)
    step!(env, x, u)
    simulate!(env; record=true)
    visualize(env)
    @test true
end
