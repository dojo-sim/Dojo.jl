using Dojo
using FiniteDiff
using JLD2
using LinearAlgebra
using Plots
using PyCall

################################################################################
# build NeRf
################################################################################
@pyinclude(joinpath(OSF_PATH, "extract_density_julia_cpu.py"))
nerf_object = py"generate_test_nerf"()

# Nerf data
deps_folder = joinpath(module_dir(), "environments/bunny/deps")
outer_mesh = jldopen(joinpath(deps_folder, "bunny_outer_mesh.jld2"))["mesh"]
inner_mesh = jldopen(joinpath(deps_folder, "bunny_inner_mesh.jld2"))["mesh"]

vis = Visualizer()
open(vis)

################################################################################
# Offline
################################################################################

halfspace_origin = [0,0,0.1]
normal = [0.2,0,1.0]
halfspace = HalfSpaceCollider(halfspace_origin, normal)
soft = SoftCollider(nerf_object, N=5000)

jldsave(joinpath(deps_folder, "soft_collider.jld2"), soft=soft)

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



# α = 0.2
# soft1 = deepcopy(soft)
# soft2 = deepcopy(soft)
# soft1.x = α*[0,0,-1.0]
# soft2.x = α*[0,0,+1.0]
# soft1.q = Quaternion(0,1,0,0.0)
# soft2.q = Quaternion(1,0,0,0.0)
# cross_collision(soft1, soft2)
# cross_collision(soft2, soft1)
# collision(soft1, soft2)
# collision(soft1, soft2)
