# ## Ghost
set_camera!(vis, cam_pos=[-1,1,0], zoom=1)

z_sim = get_maximal_state(storage)
timesteps = [5, 10, 15]# .+ 150

for t in timesteps
    name = Symbol("robot_$t")
    build_robot(mech, vis=vis, name=name, color= magenta_light)
    z = z_sim[t]
    set_robot(vis, mech, z, name=name)
end

z_sim[220][1] -= 0.1
z_sim[220][13 + 1] -= 0.05
set_robot(vis, mech, z_sim[220])

z_sim[104][1] -= 0.025
z_sim[104][13 + 1] -= 0.05
set_robot(vis, mech, z_sim[104])
