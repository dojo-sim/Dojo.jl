using Dojo
using JLD2
include("../ags.jl")

# ## Environment
env = make("halfcheetah", dt=0.05)
env.nx
obs = reset(env)

# ## Open visualizer
open(env.vis)

# ## Augmented random search

# # ## Training
# train_times = Float64[]
# rewards = Float64[]
# policies = Matrix{Float64}[]

# Reset environment
env = make("halfcheetah", dt=0.05)
obs = reset(env)

# Random policy
hp = HyperParameters(main_loop_size = 30, horizon = 80, n_directions = 6, b = 6, step_size = 0.02)
input_size = length(obs)
output_size = length(env.u_prev)
normalizer = Normalizer(input_size)
policy = Policy(input_size, output_size, hp)

# Train policy
train_time = @elapsed train(env, policy, normalizer, hp)

# Evaluate policy
reward = rollout_policy(policy.θ, env, normalizer, hp)
traj = display_policy(env, policy, normalizer, hp)
visualize(env, traj)
open(env.vis)





datamat = full_data_matrix(env.mechanism, attjac = true)

joint0 = env.mechanism.joints[3]
pbody0 = get_body(env.mechanism, joint0.parent_id)
cbody0 = get_body(env.mechanism, joint0.child_ids[1])
constraint_jacobian_parent(joint0.constraints[1], pbody0, cbody0, joint0.variables[2], env.mechanism.timestep)
constraint_jacobian_parent(joint0.constraints[2], pbody0, cbody0, joint0.variables[2], env.mechanism.timestep)


constraint_jacobian_child(joint0.constraints[1], pbody0, cbody0, joint0.variables[2], env.mechanism.timestep)
constraint_jacobian_child(joint0.constraints[2], pbody0, cbody0, joint0.variables[2], env.mechanism.timestep)
length(joint0.constraints[2])




# Cache
push!(train_times, train_time)
push!(rewards, reward)
push!(policies, policy.θ)

# N = 5
# for i = 1:N
#     # Reset environment
#     env = make("halfcheetah", dt=0.05)
#     obs = reset(env)
#
#     # Random policy
#     hp = HyperParameters(main_loop_size = 30, horizon = 80, n_directions = 6, b = 6, step_size = 0.02)
#     input_size = length(obs)
#     output_size = length(env.u_prev)
#     normalizer = Normalizer(input_size)
#     policy = Policy(input_size, output_size, hp)
#
#     # Train policy
#     train_time = @elapsed train(env, policy, normalizer, hp)
#
#     # Evaluate policy
#     reward = rollout_policy(policy.θ, env, normalizer, hp)
#
#     # Cache
#     push!(train_times, train_time)
#     push!(rewards, reward)
#     push!(policies, policy.θ)
# end

# @save joinpath(@__DIR__, "halfcheetah_rl.jld2") train_times rewards policies
@load joinpath(@__DIR__, "halfcheetah_rl.jld2") train_times rewards policies

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

# ## Visualizer policy
open(env.vis)

hp = HyperParameters(main_loop_size = 30, horizon = 80, n_directions = 6, b = 6, step_size = 0.02)
input_size = length(obs)
output_size = length(env.u_prev)
normalizer = Normalizer(input_size)

traj = display_policy(env,
    # policy,
    Policy(hp, θ),
    normalizer, hp)

for t = 1:length(traj)
    traj[t][2] += 3.25
end

visualize(env, traj)

vals = [1, 2, 3, 4]
findmax(vals, 2)
max_idx = sortperm(vals, lt=Base.isgreater)
set_camera!(env.vis, cam_pos=[0,-52,0], zoom=15)




set_floor!(env.vis, z=0.03)

# ## Animation
z = [minimal_to_maximal(env.mechanism, x) for x in traj]
z = [[z[1] for t = 1:40]..., z..., [z[end] for t = 1:40]...]
T = length(z)
anim = MeshCat.Animation(convert(Int, floor(1.0 / env.mechanism.timestep)))
build_robot(env.vis, env.mechanism, color=magenta)
for t = 1:T
    MeshCat.atframe(anim, t) do
        set_robot(env.vis, env.mechanism, z[t])
    end
end
MeshCat.setanimation!(env.vis, anim)

# ## Ghost
env = make("halfcheetah", dt=0.05)
open(env.vis)
setvisible!(env.vis[:robot], false)
timesteps = [1, 50, 60, 70, 80, 90, 100, 108, T]
for t in timesteps
    name = Symbol("robot_$t")
    color = (t == T ? magenta : magenta_light)
    build_robot(env.vis, env.mechanism, color=color, name=name)
    set_robot(env.vis, env.mechanism, z[t], name=name)
end

# ## Save/Load policy
# θ = policy.θ
# @save joinpath(@__DIR__, "halfcheetah_policy.jld2") θ
@load joinpath(@__DIR__, "halfcheetah_policy.jld2") θ


## test random policy
env = make("halfcheetah", mode=:min, g=-9.81, dt=0.05)
initialize!(env.mechanism, :halfcheetah)
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
