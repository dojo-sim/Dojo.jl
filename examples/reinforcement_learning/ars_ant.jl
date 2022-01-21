using Dojo
using JLD2 
include("ars.jl")

# ## Ant
env = make("ant", mode=:min, g=-9.81, dt=0.05, damper=10.0, spring=1.0, cf = 0.5,
    contact=true, contact_body=true)
obs = reset(env)
initializeant!(env.mechanism, pos = [1.3,0,0], rot = [0,0,0.])
env.x .= getMinState(env.mechanism)
render(env)

# ## Open visualizer
open(env.vis)

# ## Set up policy
hp = HyperParameters(main_loop_size=30, horizon=150, n_directions=6, b=6, step_size=0.02)
input_size = length(obs)
output_size = length(env.u_prev)

# ## Training 
train_times = Float64[] 
rewards = Float64[]
policies = Matrix{Float64}[]
N = 1 
for i = 1:N
    # Random policy
    normalizer = Normalizer(input_size)
    policy = Policy(input_size, output_size, hp)

    # Train policy
    train_time = @elapsed train(env, policy, normalizer, hp)

    # Evaluate policy
    reward = rollout_policy(policy.θ, env, normalizer, hp)

    # Cache 
    # push!(train_times, train_time) 
    push!(rewards, reward) 
    push!(policies, policy.θ)
end


# ## Visualize policy
# traj = display_random_policy(env, hp)
normalizer = Normalizer(input_size)
policy = Policy(hp, policies[1])

# Train policy
train_time = @elapsed train(env, policy, normalizer, hp)

traj = display_policy(env, policy, normalizer, hp)
# visualize(env, traj)

for t = 1:length(traj)
    traj[t][1] -= 1.5 
end
z = [min2max(env.mechanism, x) for x in traj]
z = [[z[1] for t = 1:40]..., z..., [z[end] for t = 1:40]...]
T = length(z) 

anim = MeshCat.Animation(convert(Int, floor(1.0 / env.mechanism.Δt)))
build_robot(env.vis, env.mechanism, color=cyan)
for t = 1:T
    MeshCat.atframe(anim, t) do
        set_robot(env.vis, env.mechanism, z[t])
    end
end
MeshCat.setanimation!(env.vis, anim)

MeshCat.settransform!(env.vis["/Cameras/default"],
        MeshCat.compose(MeshCat.Translation(0.0, 0.0, 3.0), MeshCat.LinearMap(Rotations.RotZ(-π / 2.0))))
setprop!(env.vis["/Cameras/default/rotated/<object>"], "zoom", 0.75)

# ## Ghost 
env = make("ant", mode=:min, g=-9.81, dt=0.05, damper=10.0, spring=1.0, cf = 0.5,
    contact=true, contact_body=true)
open(env.vis)
setvisible!(env.vis[:robot], false)
timesteps = [1, 70, 110, 150, T] 
for t in timesteps
    name = Symbol("robot_$t")
    color = (t == T ? cyan : cyan_light)
    build_robot(env.vis, env.mechanism, color=color, name=name)
    set_robot(env.vis, env.mechanism, z[t], name=name)
end

# ## Save/Load policy
θ = policy.θ
# @save joinpath(@__DIR__, "ant_policy.jld2") θ
# @load joinpath(@__DIR__, "ant_policy.jld2") θ

# ## test random policy
env = make("ant", mode=:min, g=-9.81, dt=0.05, damper=10.0, spring=1.0, cf = 0.5,
    contact=true, contact_body=true)
# initialize!(env.mechanism, :ant)
open(env.vis)
# storage = simulate!(env.mechanism, 1.0, record=true, verbose=false)
# visualize(env.mechanism, storage, vis=env.vis, show_contact=true)

reset(env)
render(env)
x0 = getMinState(env.mechanism)

for i = 1:100
    u = rand(Uniform(-1.0, 1.0), env.nu)
    x0, r, _ = step(env, x0, u)
    @show r
    render(env)
end