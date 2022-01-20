using Dojo
using JLD2 
include("ars.jl")

# ## Environment 
env = make("halfcheetah", dt=0.05)
obs = reset(env)

# ## Open visualizer
open(env.vis)

# ## Augmented random search
hp = HyperParameters(main_loop_size = 100, horizon = 80, n_directions = 6, b = 6, step_size = 0.02)
input_size = length(obs)
output_size = length(env.u_prev)
policy = Policy(input_size, output_size, hp)
normalizer = Normalizer(input_size)

# ## Train policy
train(env, policy, normalizer, hp)

# ## Visualizer policy
open(env.vis)

traj = display_policy(env, policy, normalizer, hp)
for t = 1:length(traj) 
    traj[t][2] += 3.25
end
visualize(env, traj)
MeshCat.settransform!(env.vis["/Cameras/default"],
        MeshCat.compose(MeshCat.LinearMap(Rotations.RotZ(-π / 2.0)), MeshCat.Translation(50.0, 0.0, -1.0)))
setprop!(env.vis["/Cameras/default/rotated/<object>"], "zoom", 15.0)


# ## Animation
z = [min2max(env.mechanism, x) for x in traj]
z = [[z[1] for t = 1:40]..., z..., [z[end] for t = 1:40]...]
T = length(z) 
anim = MeshCat.Animation(convert(Int, floor(1.0 / env.mechanism.Δt)))
build_robot(env.vis, env.mechanism, color=RGBA(1.0, 153.0 / 255.0, 51.0 / 255.0, 1.0))
for t = 1:T
    MeshCat.atframe(anim, t) do
        set_robot(env.vis, env.mechanism, z[t])
    end
end
MeshCat.setanimation!(env.vis, anim)

# ## Ghost 
env = make("halfcheetah", vis=vis, dt=0.05)
vis = Visualizer()
open(env.vis)
setvisible!(env.vis[:robot], false)
timesteps = [1, 40, 60, 70, 80, 85, 90, 95, T] 
for t in timesteps
    name = Symbol("robot_$t")
    color = (t == T ? orange : orange_light)
    build_robot(env.vis, env.mechanism, color=color, name=name)
    set_robot(env.vis, env.mechanism, z[t], name=name)
end

# ## Save/Load policy
# @save joinpath(@__DIR__, "halfcheetah_policy.jld2") θ
@load joinpath(@__DIR__, "halfcheetah_policy.jld2") θ


