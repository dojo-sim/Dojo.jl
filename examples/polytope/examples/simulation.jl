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

sphere_origin = [0,0,0.0]
radius = 0.2
sphere = SphereCollider(sphere_origin, radius)


# soft = SoftCollider(nerf_object, mesh, N=5000)
soft = jldopen(joinpath(example_dir(), "results", "soft.jld2"))["soft"]

@elapsed collision(halfspace, soft)
@elapsed collision(sphere, soft)
@elapsed collision(soft, soft)

################################################################################
# Simulation
################################################################################
set_light!(vis)
set_floor!(vis, color=RGBA(0.4,0.4,0.4,0.8), z=0.5, origin=halfspace.origin, normal=halfspace.normal)
set_background!(vis)

timestep = 0.01
H = Int(floor(0.6/timestep))
mass = soft.mass
inertia = soft.inertia
gravity = [0,0,-9.81]

x2 = [0,0,1.3]
v15 = 0*[0,2,3.0]
q2 = rand(4)
q2 /= norm(q2)
q2 = Quaternion(q2...)
q2 = Quaternion(1,0,0,0.0)
ϕ15 = 0*[3.5,0,0.0]
z0 = [x2; v15; Dojo.vector(q2); ϕ15]

ψ, impact_normal, barycenter = collision(sphere, soft)


soft.options = ColliderOptions19(
    impact_damper=1e7,
    impact_spring=3e6,
    sliding_drag=0.0,
    sliding_friction=0.30,
    rolling_drag=0.05,
    rolling_friction=0.01,
)
# @time Z0, Φ0 = implicit_simulation(halfspace, soft, timestep, H, z0)
# @time Z0, Φ0 = implicit_simulation(sphere, soft, timestep, H, z0)
soft1 = deepcopy(soft)
soft1.x = [0,0,0.0]
soft1.q = Quaternion(0,1,0,0.,true)
@time Z0, Φ0 = implicit_simulation(soft1, soft, timestep, H, z0)
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

collision(soft1, soft)

build_mesh!(tight_mesh, vis, name=:tight_mesh, color=RGBA(1,1,1,1.0))
build_mesh!(mesh, vis, name=:mesh, color=RGBA(1,1,1,0.3))
build_mesh!(tight_mesh, vis, name=:tight_mesh_1, color=RGBA(0,0,0,1.0))
build_mesh!(mesh, vis, name=:mesh_1, color=RGBA(1,1,1,0.3))
build_sphere!(sphere, vis, name=:sphere, color=RGBA(0.4,0.4,0.4,1.0))
animation = MeshCat.Animation(Int(floor(1/timestep)))
for i = 1:H
    atframe(animation, i) do
        set_robot!(Z0[i][1:3], Quaternion(Z0[i][7:10]...), vis, name=:mesh)
        set_robot!(Z0[i][1:3], Quaternion(Z0[i][7:10]...), vis, name=:tight_mesh)
        set_robot!(soft1.x, soft1.q, vis, name=:mesh_1)
        set_robot!(soft1.x, soft1.q, vis, name=:tight_mesh_1)
    end
end
setanimation!(vis, animation)
# render_static(vis)
# open(joinpath(results_dir, "bunny_hull.html"), "w") do file
#     write(file, static_html(vis))
# end
Dojo.convert_frames_to_video_and_gif("bunny_on_top_of_bunny")
