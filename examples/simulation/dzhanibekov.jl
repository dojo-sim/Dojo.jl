################################################################################
# Dzhanibekov effect
################################################################################
using Dojo
using MeshCat

# ## Visualizers
vis = Visualizer()
open(vis)

# ## Simulation
timestep = 0.01
gravity = 0.0
mech = getdzhanibekov(Δt=timestep, g=gravity);
initializedzhanibekov!(mech, ω=[15.0; 0.01; 0.0])
storage = simulate!(mech, 4.65, record=true, verbose=false)

# ## Simulation
visualize(mech, storage, vis=vis)

# ## Ghost
set_camera!(vis, cam_pos=[-1,1,0], zoom=1)

z_sim = get_max_state(storage)
timesteps = [5, 10, 15]# .+ 150

for t in timesteps
    name = Symbol("robot_$t")
    build_robot(vis, mech, name=name, color= magenta_light)
    z = z_sim[t]
    set_robot(vis, mech, z, name=name)
end

z_sim[220][1] -= 0.1
z_sim[220][13 + 1] -= 0.05
set_robot(vis, mech, z_sim[220])

z_sim[104][1] -= 0.025
z_sim[104][13 + 1] -= 0.05
set_robot(vis, mech, z_sim[104])
