using Dojo
using JLD2 
include("ars.jl")

# ## Ant
env = make("ant", vis = vis, mode=:min, g=-9.81, dt=0.05, damper=50.0, spring=30.0, cf = 0.5,
    contact=true, contact_body=true)
obs = reset(env)
initializeant!(env.mechanism, pos = [1.3,0,0], rot = [0,0,0.])
env.x .= getMinState(env.mechanism)
render(env)

# ## Open visualizer
open(env.vis)

# ## Set up policy
hp = HyperParameters(main_loop_size=100, horizon=150, n_directions=6, b=6, step_size=0.02)
input_size = length(obs)
output_size = length(env.u_prev)
policy = Policy(input_size, output_size, hp)
normalizer = Normalizer(input_size)

# ## Train policy
3
train(env, policy, normalizer, hp)

# ## Visualize policy
# traj = display_random_policy(env, hp)
traj = display_policy(env, policy, normalizer, hp)
# visualize(env, traj)

for t = 1:length(traj)
    traj[t][1] -= 1.5 
end
z = [min2max(env.mechanism, x) for x in traj]
z = [[z[1] for t = 1:40]..., z..., [z[end] for t = 1:40]...]
T = length(z) 

anim = MeshCat.Animation(convert(Int, floor(1.0 / env.mechanism.Δt)))
build_robot(env.vis, env.mechanism, color=RGBA(51.0 / 255.0, 1.0, 1.0, 1.0))
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
env = make("ant", vis = vis, mode=:min, g=-9.81, dt=0.05, damper=50.0, spring=30.0, cf = 0.5,
    contact=true, contact_body=true)
vis = Visualizer()
open(env.vis)
setvisible!(env.vis[:robot], false)
timesteps = [1, 70, 110, 150, T] 
for t in timesteps
    name = Symbol("robot_$t")
    color = (t == T ? RGBA(51.0 / 255.0, 1.0, 1.0, 1.0) : RGBA(51.0 / 255.0, 1.0, 1.0, 0.25))
    build_robot(env.vis, env.mechanism, color=color, name=name)
    set_robot(env.vis, env.mechanism, z[t], name=name)
end

# ## Save/Load policy
θ = policy.θ
@save joinpath(@__DIR__, "ant_policy.jld2") θ
@load joinpath(@__DIR__, "ant_policy.jld2") θ
