# Open visualizer
env = make("ant", mode=:min, g=-9.81, dt=0.05, damper=25.0, spring=100.0, contact=true, contact_body=true);
obs = reset(env)
render(env)
open(env.vis)
hp = HyperParameters(env_name="ant", main_loop_size = 30, horizon = 100, n_directions = 8, b = 8, step_size = 0.02)
input_size = length(obs)
output_size = length(env.u_prev)
policy = Policy(input_size, output_size, hp)
normalizer = Normalizer(input_size)

train(env, policy, normalizer, hp)

traj = display_policy(env, policy, normalizer, hp)
visualize(env, traj)
