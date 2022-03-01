using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# ## Setup
using Dojo

# ## Visualizer
vis= Visualizer()
open(vis)
set_camera!(vis,
    cam_pos=[-0.01, 0.0, 90.0],
    zoom=30)
set_floor!(vis,
    x=20.0,
    y=20.0,
    color=RGBA(1.0, 1.0, 1.0, 1.0))
set_light!(vis,
    ambient=0.65,
    direction="Positive")

################################################################################
# Nonlinear Friction Cone vs. Linearized Friction Cone
################################################################################
timestep=0.01
gravity=-9.81
friction_coefficient = 0.25
x0 = [-1.5, -0.50, 0.25]
v0 = [4, 0.80, 0.0]
ω0 = [0.0, 0.0, 0.0]
opts = SolverOptions(rtol=1.0e-6, btol=1.0e-6)

# ## Linear cone
color_lc = orange;
mech_lc = get_mechanism(:box,
    timestep=timestep,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    contact_type=:linear,
    mode=:box,
    color=color_lc);

# ## Simulate
initialize!(mech_lc, :box,
    x=x0,
    q=one(UnitQuaternion),
    v=v0,
    ω=ω0)

storage_lc = simulate!(mech_lc, 4.0,
    record=true,
    opts=opts)

# ## Visualize
via, anim = visualize(mech_lc, storage_lc,
    vis=vis,
    name=:lc)

line_mat_lc = LineBasicMaterial(color=color_lc, linewidth=10.0)
points_lc = Vector{Point{3,Float64}}()
for (i, xt) in enumerate(storage_lc.x[1])
    k = xt
    push!(points_lc, Point(k...))
end
setobject!(vis[:path_lc], MeshCat.Line(points_lc, line_mat_lc))

# ## Nonlinear cone
color_nc = cyan;
mech_nc = get_mechanism(:box,
    timestep=timestep,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    contact_type=:nonlinear,
    mode=:box,
    color=color_nc);

# ## Simulate
initialize!(mech_nc, :box,
    x=x0,
    q=one(UnitQuaternion),
    v=v0,
    ω=ω0)

storage_nc = simulate!(mech_nc, 4.0,
    record=true,
    opts=opts)

# ## Visualize
line_mat_nc = LineBasicMaterial(color=color_nc, linewidth=25.0)
points_nc = Vector{Point{3,Float64}}()
for (i, xt) in enumerate(storage_nc.x[1])
    k = xt
    @show k[3]
    push!(points_nc, Point(k...))
end
setobject!(vis[:path_nc], MeshCat.Line(points_nc, line_mat_nc))

vis, anim = visualize(mech_nc, storage_nc,
    vis=vis,
    name=:nc,
    animation=anim)

# ## MuJoCo pyramidal cone
color_mjlc = magenta;

mech_mjlc = get_mechanism(:box,
    timestep=timestep,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    contact_type=:linear,
    mode=:box, color=color_mjlc);

# ## Load
initialize!(mech_mjlc, :box,
    x=x0,
    q=one(UnitQuaternion),
    v=v0,
    ω=ω0)
file = jldopen(joinpath(@__DIR__, "../mujoco_benchmark/results/cone_compare_pyramidal.jld2"))
storage_mjlc = generate_storage(mech_mjlc, [get_maximal_state(mech_mjlc), file["ztraj"]...])

# ## Visualize
vis, anim = visualize(mech_mjlc, storage_mjlc,
    vis=vis,
    name=:mjlc,
    animation=anim)

line_mat_mjlc = LineBasicMaterial(color=color_mjlc, linewidth=25.0)
points_mjlc = Vector{Point{3,Float64}}()
for (i, xt) in enumerate(storage_mjlc.x[1])
    k = xt
    push!(points_mjlc, Point(k[1], k[2], k[3]))
end
setobject!(vis[:path_mjlc], MeshCat.Line(points_mjlc, line_mat_mjlc))


# ## MuJoCo elliptic cone
color_mjnc = RGBA(0,0,0);

mech_mjnc = get_mechanism(:box,
    timestep=timestep,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    contact_type=:linear,
    mode=:box, color=color_mjnc);

# ## Load
initialize!(mech_mjnc, :box,
    x=x0,
    q=one(UnitQuaternion),
    v=v0,
    ω=ω0)
file = jldopen(joinpath(@__DIR__, "../mujoco_benchmark/results/cone_compare_elliptic.jld2"))
storage_mjnc = generate_storage(mech_mjnc, [get_maximal_state(mech_mjnc), file["ztraj"]...])

# ## Visualize
vis, anim = visualize(mech_mjnc, storage_mjnc,
    vis=vis,
    name=:mjnc,
    animation=anim)

line_mat_mjnc = LineBasicMaterial(color=color_mjnc, linewidth=25.0)
points_mjnc = Vector{Point{3,Float64}}()
for (i, xt) in enumerate(storage_mjnc.x[1])
    k = xt
    push!(points_mjnc, Point(k[1], k[2], k[3]))
end
setobject!(vis[:path_mjnc], MeshCat.Line(points_mjnc, line_mat_mjnc))

settransform!(vis[:lc], MeshCat.Translation(0,+0.04,0))
settransform!(vis[:path_lc], MeshCat.Translation(0,+0.04,0))
settransform!(vis[:nc], MeshCat.Translation(0,-0.04,0))
settransform!(vis[:path_nc], MeshCat.Translation(0,-0.04,0))
settransform!(vis[:mjlc], MeshCat.Translation(0,+0.00,0))
# settransform!(vis[:mjnc], MeshCat.Translation(0,-0.08,-0.02))
settransform!(vis[:mjnc], MeshCat.Translation(0,-0.08,0))
settransform!(vis[:path_mjnc], MeshCat.Translation(0,-0.08,0))
