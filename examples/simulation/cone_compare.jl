# ## visualizer
vis = Visualizer()
open(vis)

MeshCat.settransform!(vis["/Cameras/default"],
        MeshCat.compose(MeshCat.Translation(0.0, 0.0, 50.0), MeshCat.LinearMap(Rotations.RotY(-pi / 2.5))))
setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 10)

################################################################################
# Nonlinear Friction Cone vs. Linearized Friction Cone
################################################################################
timestep = 0.01
gravity = -9.81
friction_coefficient = 1.0

# ## linear cone
color_lc = RGBA(1.0, 153.0 / 255.0, 51.0 / 255.0, 1.0);
mech_lc = getmechanism(:box, Δt=timestep, g=gravity, cf=friction_coefficient, conetype=:linear, mode=:box, color=color_lc);
initialize!(mech_lc, :box, x=[-3.0, 0.0, 1.0], q=one(UnitQuaternion), 
    v=[10.0, 2.0, 0.0], ω=[0.0, 0.0, 0.0])
storage_lc = simulate!(mech_lc, 5.0, record=true, verbose=false)
visualize(mech_lc, storage_lc, vis=vis, name=:lc)

line_mat_lc = LineBasicMaterial(color=color_lc, linewidth=10.0)
points_lc = Vector{Point{3,Float64}}()
for (i, xt) in enumerate(storage_lc.x[1])
    k = xt
    push!(points_lc, Point(k[1], k[2], k[3]))
end
setobject!(vis[:path_lc], MeshCat.Line(points_lc, line_mat_lc))

# ## nonlinear cone
color_nc = RGBA(51.0 / 255.0, 1.0, 1.0, 1.0);
mech_nc = getmechanism(:box, Δt=timestep, g=gravity, cf=friction_coefficient, conetype=:soc, mode=:box, color=color_nc);
initialize!(mech_nc, :box, x=[-3.0, 0.0, 1.0], q=one(UnitQuaternion), 
    v=[10.0, 2.0, 0.0], ω=[0.0, 0.0, 0.0])
storage_nc = simulate!(mech_nc, 5.0, record=true, verbose=false)
visualize(mech_nc, storage_nc, vis=vis, name=:nc)

line_mat_nc = LineBasicMaterial(color=color_nc, linewidth=25.0)
points_nc = Vector{Point{3,Float64}}()
for (i, xt) in enumerate(storage_nc.x[1])
    k = xt
    push!(points_nc, Point(k[1], k[2], k[3] + 0.05))
end
setobject!(vis[:path_nc], MeshCat.Line(points_nc, line_mat_nc))

# ## joint animation 
max_lc = getMaxState(storage_lc)
max_nc = getMaxState(storage_nc)
max_lc = [[max_lc[1] for t = 1:100]..., max_lc...]
max_nc = [[max_nc[1] for t = 1:100]..., max_nc...]
T = length(max_nc) 
anim = MeshCat.Animation(convert(Int, floor(1.0 / timestep)))
for t = 1:T
    MeshCat.atframe(anim, t) do
        set_robot(vis, mech_lc, max_lc[t]; name=:lc)
        set_robot(vis, mech_nc, max_nc[t]; name=:nc)
    end
end
MeshCat.setanimation!(vis, anim)
