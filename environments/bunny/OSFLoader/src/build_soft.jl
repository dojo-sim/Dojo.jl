using Dojo
using OSFLoader
using PyCall

# run Dojo's activate.jl
# run using Dojo
# run OSFLoader's activate.jl
# run using OSFLoader, PyCall
# rerun Dojo's activate.jl
# this script needs to be executed in the Dojo module and Dojo Env after loading both Dojo and the OSFLOader modules

################################################################################
# build NeRf
################################################################################
nerf_object = get_nerf_object()

# Nerf data
deps_folder = joinpath(module_dir(), "environments/bunny/deps")
outer_mesh = jldopen(joinpath(deps_folder, "bunny_outer_mesh.jld2"))["mesh"]
inner_mesh = jldopen(joinpath(deps_folder, "bunny_inner_mesh.jld2"))["mesh"]

vis = Visualizer()
open(vis)

################################################################################
# Offline
################################################################################

soft = SoftCollider(nerf_object, N=5000)
jldsave(joinpath(deps_folder, "soft_collider.jld2"), soft=soft)
