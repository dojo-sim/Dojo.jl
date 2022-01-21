using Dojo

# ## Open visualizer
vis = Visualizer()
open(vis)

# ## Mechanism
mech = getmechanism(:atlas, Δt=0.01, g=-9.81, cf=0.5, damper=100.0, spring=1.0, contact=true)
@show length(mech.bodies) * 13
@show controldim(mech)

# ## Simulate
initializeatlasstance!(mech, tran=[0,0,0.5], rot=[0.0,0.0, 0.0])
storage = simulate!(mech, 2.25, record=true, opts=InteriorPointOptions(rtol=1.0e-6, btol=1e-6))

# ## Visualize
visualize(mech, storage, vis=vis)

# ## Animation
z = getMaxState(storage)
z = [[z[1] for t = 1:100]..., z..., [z[end] for t = 1:100]...]
T = length(z)
anim = MeshCat.Animation(convert(Int, floor(1.0 / 0.01)))
build_robot(vis, mech)
for t = 1:T
    MeshCat.atframe(anim, t) do
        set_robot(vis, mech, z[t])
    end
end
MeshCat.setanimation!(vis, anim)
MeshCat.settransform!(vis["/Cameras/default"],
        MeshCat.compose(MeshCat.Translation(0.0, 0.0, -1.002), MeshCat.LinearMap(Rotations.RotY(0.0))))
setprop!(vis["/Cameras/default/rotated/<object>"], "zoom",0.840)

set_floor!(vis, x=0.2, y=10, z=0.02, color=RGBA(0,0,0,1))
