using LinearAlgebra
using Plots
using FiniteDiff
using JLD2
using Dojo
using PyCall

################################################################################
# build NeRf
################################################################################
@pyinclude(joinpath(OSF_PATH, "extract_density_julia.py"))
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

halfspace_origin = [0,0,0.1]
normal = [0.2,0,1.0]
halfspace = HalfSpaceCollider(halfspace_origin, normal)
soft = SoftCollider(nerf_object, mesh, N=5000)

# jldsave(joinpath(example_dir(), "results", "soft.jld2"), soft=soft)

@time collision(halfspace, soft)

X = -1.0:0.01:2.0
phi = zeros(length(X))
∇phi = zeros(length(X),3)
for (i,x) in enumerate(X)
    soft.x = [0,0,x]
    phi[i], ∇phi[i,:], Nphi[i,:] = collision(halfspace, soft)
end
plot(X, phi)
plot!(X, ∇phi[:,1])
plot!(X, ∇phi[:,2])
plot!(X, ∇phi[:,3])

set_light!(vis)
set_background!(vis)
set_floor!(vis, alt=-0.5)
build_collider!(soft, vis, visualize_particle=false)
