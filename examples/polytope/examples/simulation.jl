################################################################################
# Translation only
################################################################################
using LinearAlgebra
using Plots
using FiniteDiff
using Dojo
using PyCall

################################################################################
# build NeRf
################################################################################
@pyinclude(joinpath(OSF_PATH, "extract_density_julia_cpu.py"))
nerf_object = py"generate_test_nerf"()

# Nerf data
results_dir = joinpath(example_dir(), "results")
mesh = jldopen(joinpath(results_dir, "bunny_mesh.jld2"))["mesh"]
tight_mesh = jldopen(joinpath(results_dir, "bunny_tight_mesh.jld2"))["mesh"]

vis = Visualizer()
open(vis)


################################################################################
# Offline
################################################################################
halfspace_origin = [0,0,-0.5]
normal = [0,0,1.0]
halfspace = HalfSpaceCollider(halfspace_origin, normal)

sphere_origin = [0,0,-5.0]
radius = 5.2
sphere = SphereCollider(sphere_origin, radius)


# soft0 = SoftCollider(nerf_object, mesh, N=5000)
# jldsave(joinpath(example_dir(), "results", "soft.jld2"), soft=soft0)
soft0 = jldopen(joinpath(example_dir(), "results", "soft.jld2"))["soft"]
soft1 = jldopen(joinpath(example_dir(), "results", "soft.jld2"))["soft"]
soft1.x = [0,0,-0.13]
soft1.q = Quaternion(1,0,0,0.,true)


@elapsed collision(halfspace, soft0)
@elapsed collision(sphere, soft0)
@elapsed collision(soft0, soft0)

################################################################################
# Simulation
################################################################################
set_light!(vis)
set_floor!(vis, color=RGBA(0.4,0.4,0.4,0.8), z=0.5, origin=halfspace.origin, normal=halfspace.normal)
set_background!(vis)

timestep = 0.01
H = Int(floor(1.5/timestep))
gravity = [0,0,-9.81]

x2 = [0,-1.5,0.5]
v15 = [0,5,3.0]
q2 = rand(4)
# q2 = [1,0.05,0,0]
q2 = Quaternion(q2 ./ norm(q2)...,true)
ϕ15 = [3.5,0,0.0]
z0 = [x2; v15; Dojo.vector(q2); ϕ15]

ψ, barycenter, contact_normal = collision(sphere, soft0)
ψ, barycenter, contact_normal = collision(halfspace, soft0)
ψ, barycenter, contact_normal = collision(soft1, soft0)

soft0.options = ColliderOptions110( # nerf nerf
    impact_damper=1e7,
    impact_spring=1e7,
    sliding_drag=0.0,
    sliding_friction=0.1,
    rolling_drag=0.05,
    rolling_friction=0.01,
)
# soft0.options = ColliderOptions110( # halfspace
#     impact_damper=3e4,
#     impact_spring=1e4,
#     sliding_drag=0.0,
#     sliding_friction=0.05,
#     rolling_drag=0.02,
#     rolling_friction=0.02,
# )
# @time Z0, Φ0 = implicit_simulation(halfspace, soft0, timestep, H, z0, gravity)
# @time Z0, Φ0 = implicit_simulation(sphere, soft0, timestep, H, z0, gravity)
@time Z0, Φ0 = implicit_simulation(soft1, soft0, timestep, H, z0, gravity)

plt = plot(layout=(3,1), xlabel="time")
plot!(plt[1,1], [i*timestep for i=0:H], [z[1] for z in Z0], linewidth=3.0, label="x", color=:blue)
plot!(plt[1,1], [i*timestep for i=0:H], [z[2] for z in Z0], linewidth=3.0, label="y", color=:blue)
plot!(plt[1,1], [i*timestep for i=0:H], [z[3] for z in Z0], linewidth=3.0, label="z", color=:blue)
plot!(plt[2,1], [i*timestep for i=0:H], [z[7] for
 z in Z0], linewidth=3.0, label="q1", color=:green)
plot!(plt[2,1], [i*timestep for i=0:H], [z[8] for z in Z0], linewidth=3.0, label="q2", color=:green)
plot!(plt[2,1], [i*timestep for i=0:H], [z[9] for z in Z0], linewidth=3.0, label="q3", color=:green)
plot!(plt[2,1], [i*timestep for i=0:H], [z[10] for z in Z0], linewidth=3.0, label="q4", color=:green)
plot!(plt[3,1], [i*timestep for i=1:H], Φ0, linewidth=3.0, label="ψ", color=:red)


build_mesh!(tight_mesh, vis, name=:tight_mesh, color=RGBA(1,1,1,1.0))
build_mesh!(mesh, vis, name=:mesh, color=RGBA(1,1,1,0.3))
build_mesh!(tight_mesh, vis, name=:tight_mesh_1, color=RGBA(0.2,0.2,0.2,1.0))
build_mesh!(mesh, vis, name=:mesh_1, color=RGBA(1,1,1,0.3))
build_sphere!(sphere, vis, name=:sphere, color=RGBA(0.4,0.4,0.4,1.0))
animation = MeshCat.Animation(Int(floor(1/timestep)))
[z[3] for z in Z0]
[z[1] for z in Z0]
for i = 1:H
    atframe(animation, i) do
        set_robot!(Z0[i][1:3]-vector_rotate(soft0.center_of_mass, Quaternion(Z0[i][7:10]...)), Quaternion(Z0[i][7:10]...), vis, name=:mesh)
        set_robot!(Z0[i][1:3]-vector_rotate(soft0.center_of_mass, Quaternion(Z0[i][7:10]...)), Quaternion(Z0[i][7:10]...), vis, name=:tight_mesh)
        set_robot!(soft1.x - soft1.center_of_mass, soft1.q, vis, name=:mesh_1)
        set_robot!(soft1.x - soft1.center_of_mass, soft1.q, vis, name=:tight_mesh_1)
    end
end
setanimation!(vis, animation)
# render_static(vis)
# open(joinpath(results_dir, "bunny_hull.html"), "w") do file
#     write(file, static_html(vis))
# end
# Dojo.convert_frames_to_video_and_gif("bunny_vs_bunny")