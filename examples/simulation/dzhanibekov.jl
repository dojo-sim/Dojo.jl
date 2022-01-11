################################################################################
# Dzhanibekov effect
################################################################################
using Dojo 

# ## Visualizers
vis = Visualizer()
open(vis)

# ## Simulation
timestep = 0.01
gravity = 0.0
mech = getdzhanibekov(Δt=timestep, g=gravity);
initializedzhanibekov!(mech, ω=[15.0; 0.01; 0.0])
storage = simulate!(mech, 3.75, record=true, verbose=false)

# ## Simulation
visualize(mech, storage, vis=vis)

# ## Ghost
MeshCat.settransform!(vis["/Cameras/default"],
        MeshCat.compose(MeshCat.Translation(2.0, 1.0, -1.0), MeshCat.LinearMap(Rotations.RotZ(1.0 * pi))))
setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 1)
z_sim = getMaxState(storage)
timesteps = [1, 5, 10, 15, 20] .+ 150

for t in timesteps 
    name = Symbol("robot_$t")
    build_robot(vis, mech, name=name, color=(t == timesteps[end] ? nothing : RGBA(1.0, 0.0, 0.0, 0.25)))
    z = z_sim[t] 
    set_robot(vis, mech, z, name=name)
end
