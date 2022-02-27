# ## Animation
z = get_maximal_state(storage)
z = [[z[1] for t = 1:100]..., z..., [z[end] for t = 1:100]...]
T = length(z)
anim = MeshCat.Animation(convert(Int, floor(1.0 / 0.01)))
build_robot(mech, vis=vis)
for t = 1:T
    MeshCat.atframe(anim, t) do
        set_robot(vis, mech, z[t])
    end
end
MeshCat.setanimation!(vis, anim)
set_camera!(vis, cam_pos=[3,-0,0], zoom=0.84)

set_floor!(vis, x=0.2, y=10, z=0.02, color=RGBA(0,0,0,1))

# ## HMTL scene
render_static(vis)
open(joinpath(@__DIR__, "atlas_drop.html"), "w") do file
    write(file, static_html(vis))
end
