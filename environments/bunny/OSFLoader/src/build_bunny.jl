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
bunny_nerf = get_nerf_object(filename="bunny_trans")
halfsoap_nerf = get_nerf_object(filename="halfsoap")
halfsoap_nerf["use_lightdirs"] = true
halfsoap_nerf


# Nerf data
deps_folder = joinpath(module_dir(), "environments/bunny/deps")
outer_mesh = jldopen(joinpath(deps_folder, "bunny_outer_mesh.jld2"))["mesh"]
inner_mesh = jldopen(joinpath(deps_folder, "bunny_inner_mesh.jld2"))["mesh"]

################################################################################
# Offline
################################################################################

bunny = SoftCollider(bunny_nerf, N=5000)
jldsave(joinpath(deps_folder, "bunny_collider.jld2"), collider=bunny)
halfsoap = SoftCollider(halfsoap_nerf, N=5000)

particles = zeros(Float32, 10, 3)
masses = OSFLoader.py"density_query"(halfsoap_nerf, particles)
