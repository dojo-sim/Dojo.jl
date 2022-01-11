using Dojo

# ## Open visualizer
vis = Visualizer()
open(vis)

# ## Mechanism
mech = getmechanism(:atlas, Δt=0.01, g=-9.81, cf=0.5, damper=100.0, spring=1.0, contact=true)

# ## Simulate
initialize!(mech, :atlas, tran=[0,0,0.5], rot=[0.0,0.0, 0.0])
storage = simulate!(mech, 1.75, record=true, opts=InteriorPointOptions(rtol=1.0e-6, btol=1e-6))

# ## Visualize
visualize(mech, storage, vis=vis)

# black
z = getMaxState(storage) 
z = [[z[1] for t = 1:75]..., z...]
T = length(z) 
build_robot(vis, mech)#, color=RGBA(35 / 255, 31 / 255, 32 / 255, 1.0))
anim = MeshCat.Animation(convert(Int, floor(1.0 / 0.01)))

for t = 1:T 
    MeshCat.atframe(anim, t) do
        set_robot(vis, mech, z[t])
    end
end
MeshCat.setanimation!(vis, anim)

t = 225 #1, 100, 120
set_robot(vis, mech, z[t])

MeshCat.settransform!(vis["/Cameras/default"],
        MeshCat.compose(MeshCat.Translation(0.0, 0.0, -1.0), MeshCat.LinearMap(Rotations.RotY(0.0))))
setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 0.75)