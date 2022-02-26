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
