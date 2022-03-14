using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
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
timestep = 0.01
gravity = -9.81
friction_coefficient = 0.25
x0 = [-1.5, -0.50, 0.25]
v0 = [4, 0.80, 0.0]
ω0 = [0.0, 0.0, 0.0]
opts = SolverOptions(rtol=1.0e-6, btol=1.0e-6)

# ## Linear cone
color_lc = orange;
mech_lc = get_mechanism(:block,
    timestep=timestep,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    contact_type=:linear,
    mode=:box,
    color=color_lc);

# ## Simulate
initialize!(mech_lc, :block,
    x=x0,
    q=one(Quaternion),
    v=v0,
    ω=ω0)

storage_lc = simulate!(mech_lc, 4.0,
    record=true,
    opts=opts)

# ## Visualize
via, anim = visualize(mech_lc, storage_lc,
    vis=vis,
    name=:lc)


# ## Nonlinear cone
color_nc = cyan;
mech_nc = get_mechanism(:block,
    timestep=timestep,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    contact_type=:nonlinear,
    mode=:box,
    color=color_nc);

# ## Simulate
initialize!(mech_nc, :block,
    x=x0,
    q=one(Quaternion),
    v=v0,
    ω=ω0)

storage_nc = simulate!(mech_nc, 4.0,
    record=true,
    opts=opts)

# ## Visualize
for (i, x) in enumerate(storage_nc.x[1])
    storage_nc.x[1][i] += [0.0; 0.0; 0.1]
end

vis, anim = visualize(mech_nc, storage_nc,
    vis=vis,
    name=:nc,
    animation=anim)

# ## MuJoCo pyramidal cone
color_mjlc = magenta;

mech_mjlc = get_mechanism(:block,
    timestep=timestep,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    contact_type=:linear,
    mode=:box, color=color_mjlc);

# ## Load
initialize!(mech_mjlc, :block,
    x=x0,
    q=one(Quaternion),
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

mech_mjnc = get_mechanism(:block,
    timestep=timestep,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    contact_type=:linear,
    mode=:box, color=color_mjnc);

# ## Load
initialize!(mech_mjnc, :block,
    x=x0,
    q=one(Quaternion),
    v=v0,
    ω=ω0)
file = jldopen(joinpath(@__DIR__, "../mujoco_benchmark/results/cone_compare_elliptic.jld2"))
storage_mjnc = generate_storage(mech_mjnc, [get_maximal_state(mech_mjnc), file["ztraj"]...])

# ## Visualize
vis, anim = visualize(mech_mjnc, storage_mjnc,
    vis=vis,
    name=:mjnc,
    animation=anim)

