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
env = make("ant", mode=:min, g=-9.81, dt=0.05, damper=50.0, spring=25.0, friction_coefficient = 0.5,
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


# ## test random policy
env = make("ant", mode=:min, g=-9.81, dt=0.05, damper=25.0, spring=10.0, friction_coefficient = 0.5,
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
