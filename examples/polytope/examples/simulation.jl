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
# @pyinclude(joinpath(OSF_PATH, "extract_density_julia.py"))
# nerf_object = py"generate_test_nerf"()

# Nerf data
results_dir = joinpath(example_dir(), "results")
mesh = jldopen(joinpath(results_dir, "bunny_mesh.jld2"))["mesh"]
tight_mesh = jldopen(joinpath(results_dir, "bunny_tight_mesh.jld2"))["mesh"]

vis = Visualizer()
open(vis)


################################################################################
# Offline
################################################################################
halfspace_origin = [0,0,-0.0]
normal = [0.5,0,1.0]
halfspace = HalfSpaceCollider(halfspace_origin, normal)
# soft = SoftCollider(nerf_object, mesh, N=5000)
soft = jldopen(joinpath(example_dir(), "results", "soft.jld2"))["soft"]
@time collision(halfspace, soft)

################################################################################
# Simulation
################################################################################
set_light!(vis)
set_floor!(vis, color=RGBA(0.4,0.4,0.4,0.8), z=0.5, origin=halfspace.origin, normal=halfspace.normal)
set_background!(vis)

timestep = 0.01
H = Int(floor(4.0/timestep))
mass = soft.mass
inertia = soft.inertia
gravity = [0,0,-9.81]

x2 = [0,0,1.0]
v15 = [0,2,3.0]
q2 = rand(4)
q2 ./= norm(q2)
q2 = Quaternion(q2...)
ϕ15 = [3.5,0,0.0]
z0 = [x2; v15; Dojo.vector(q2); ϕ15]

soft.options = ColliderOptions19(
    impact_damper=100.0,
    impact_spring=100.0,
    sliding_drag=0.0,
    sliding_friction=0.6,
    rolling_drag=0.05,
    rolling_friction=0.01,
)
@time Z0, Φ0 = implicit_simulation(halfspace, soft, timestep, H, z0)
plt = plot(layout=(3,1), xlabel="time")
plot!(plt[1,1], [i*h for i=0:H], [z[1] for z in Z0], linewidth=3.0, label="x", color=:blue)
plot!(plt[1,1], [i*h for i=0:H], [z[2] for z in Z0], linewidth=3.0, label="y", color=:blue)
plot!(plt[1,1], [i*h for i=0:H], [z[3] for z in Z0], linewidth=3.0, label="z", color=:blue)
plot!(plt[2,1], [i*h for i=0:H], [z[7] for z in Z0], linewidth=3.0, label="q1", color=:green)
plot!(plt[2,1], [i*h for i=0:H], [z[8] for z in Z0], linewidth=3.0, label="q2", color=:green)
plot!(plt[2,1], [i*h for i=0:H], [z[9] for z in Z0], linewidth=3.0, label="q3", color=:green)
plot!(plt[2,1], [i*h for i=0:H], [z[10] for z in Z0], linewidth=3.0, label="q4", color=:green)
plot!(plt[3,1], [i*h for i=1:H], Φ0, linewidth=3.0, label="ψ", color=:red)



build_mesh!(tight_mesh, vis, name=:tight_mesh, color=RGBA(0.6,0.8,0.9,1.0))
build_mesh!(mesh, vis, name=:mesh, color=RGBA(1,1,1,0.3))
# build_polytope!(hull, vis, name=:hull, color=RGBA(1,0.3,0.3,0.3))
animation = MeshCat.Animation(Int(floor(1/timestep)))
for i = 1:H
    atframe(animation, i) do
        set_robot!(Z0[i][1:3], Quaternion(Z0[i][7:10]...), vis, name=:mesh)
        set_robot!(Z0[i][1:3], Quaternion(Z0[i][7:10]...), vis, name=:tight_mesh)
    end
end
setanimation!(vis, animation)
# render_static(vis)
# open(joinpath(results_dir, "bunny_hull.html"), "w") do file
#     write(file, static_html(vis))
# end
# Dojo.convert_frames_to_video_and_gif("bunny_rolling_friction_minus")
