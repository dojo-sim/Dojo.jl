using Dojo
using JLD2
using Random
include("ars.jl")

# ## Ant
env = make("ant", mode=:min, g=-9.81, dt=0.05, damper=50.0, spring=25.0, cf = 0.5,
    contact=true, contact_body=true)
obs = reset(env)
initialize_ant!(env.mechanism, pos = [1.3,0,0], rot = [0,0,0.])
env.x .= get_minimal_state(env.mechanism)
render(env)

env.nx

# ## Open visualizer
open(env.vis)

# ## Set up policy
hp = HyperParameters(main_loop_size=100, horizon=150, n_directions=6, b=6, step_size=0.02)
input_size = length(obs)
output_size = length(env.u_prev)
normalizer = Normalizer(input_size)

# ## Training
train_times = Float64[]
rewards = Float64[]
policies = Matrix{Float64}[]
N = 5
for i = 1:N
    # Reset environment
    env = make("ant", mode=:min, g=-9.81, dt=0.05, damper=50.0, spring=25.0, cf = 0.5,
        contact=true, contact_body=true)
    obs = reset(env)

    # Random policy
    Random.seed!(i)
    hp = HyperParameters(main_loop_size=100, horizon=150, n_directions=6, b=6, step_size=0.02)
    input_size = length(obs)
    output_size = length(env.u_prev)
    normalizer = Normalizer(input_size)
    policy = Policy(input_size, output_size, hp)

    # Train policy
    train_time = @elapsed train(env, policy, normalizer, hp)

    # Evaluate policy
    reward = rollout_policy(policy.θ, env, normalizer, hp)

    # Cache
    push!(train_times, train_time)
    push!(rewards, reward)
    push!(policies, policy.θ)
end

# @save joinpath(@__DIR__, "ant_rl.jld2") train_times rewards policies

# Training statistics
N_best = 3
@show rewards
max_idx = sortperm(rewards, lt=Base.isgreater)
train_time_best = (train_times[max_idx])[1:N_best]
rewards_best = (rewards[max_idx])[1:N_best]
policies_best = (policies[max_idx])[1:N_best]
@show mean(train_time_best)
@show std(train_time_best)
@show mean(rewards)
@show std(rewards)

# ## Visualize policy
# traj = display_random_policy(env, hp)
normalizer = Normalizer(input_size)
policy = Policy(hp, policies_best[1])

# Train policy
train_time = @elapsed train(env, policy, normalizer, hp)

traj = display_policy(env, policy, normalizer, hp)
# visualize(env, traj)

for t = 1:length(traj)
    traj[t][1] -= 2.0
end
z = [minimal_to_maximal(env.mechanism, x) for x in traj]
z = [[z[1] for t = 1:40]..., z..., [z[end] for t = 1:40]...]
T = length(z)

anim = MeshCat.Animation(convert(Int, floor(1.0 / env.mechanism.timestep)))
build_robot(env.vis, env.mechanism, color=cyan)
for t = 1:T
    MeshCat.atframe(anim, t) do
        set_robot(env.vis, env.mechanism, z[t])
    end
end
MeshCat.setanimation!(env.vis, anim)
set_camera!(env.vis, cam_pos=[0,0,90], zoom=20)




# ## Ghost
env = make("ant", mode=:min, g=-9.81, dt=0.05, damper=50.0, spring=25.0, cf = 0.5,
    contact=true, contact_body=true)
open(env.vis)
setvisible!(env.vis[:robot], false)
timesteps = [1, 70, 110, 130, 150, T]
for t in timesteps
    name = Symbol("robot_$t")
    color = (t == T ? cyan : cyan_light)
    build_robot(env.vis, env.mechanism, color=color, name=name)
    set_robot(env.vis, env.mechanism, z[t], name=name)
end

# ## Save/Load policy
# θ = policy.θ
# @save joinpath(@__DIR__, "ant_policy.jld2") θ
# @load joinpath(@__DIR__, "ant_policy.jld2") θ

# ## test random policy
env = make("ant", mode=:min, g=-9.81, dt=0.05, damper=25.0, spring=10.0, cf = 0.5,
    contact=true, contact_body=true)
# initialize!(env.mechanism, :ant)
open(env.vis)
# storage = simulate!(env.mechanism, 1.0, record=true, verbose=false)
# visualize(env.mechanism, storage, vis=env.vis, show_contact=true)

reset(env)
render(env)
x0 = get_minimal_state(env.mechanism)

for i = 1:100
    u = rand(Uniform(-1.0, 1.0), env.nu)
    x0, r, _ = step(env, x0, u)
    @show r
    render(env)
end
