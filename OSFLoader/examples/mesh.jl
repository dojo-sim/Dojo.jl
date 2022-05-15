using OSFLoader
using GeometryBasics
using LinearAlgebra
using JLD2
using MeshIO
using Meshing
using MeshCat
using Quaternions

vis = Visualizer()
open(vis)

results_dir = joinpath(OSFLoader.osf_loader_dir(), "assets/mesh")

################################################################################
# bunny
################################################################################
bunny_nerf = OSFLoader.get_nerf_object(filename="bunny")
mesh_low, mesh_high = nerf_mesh(bunny_nerf, sampling_density=100,
        normalizer=DensityFieldNormalizer(nerf=:bunny))
setobject!(vis[:low], mesh_low, MeshPhongMaterial(color=Colors.RGBA(1,1,1,0.5)))
setobject!(vis[:high], mesh_high, MeshPhongMaterial(color=Colors.RGBA(0,0,0,1)))

save(joinpath(results_dir, "bunny_low.obj"), mesh_low)
save(joinpath(results_dir, "bunny_high.obj"), mesh_high)

################################################################################
# bluesoap
################################################################################
bluesoap_nerf = OSFLoader.get_nerf_object(filename="bluesoap")
mesh_low, mesh_high = nerf_mesh(bluesoap_nerf, sampling_density=100,
        normalizer=DensityFieldNormalizer(nerf=:bluesoap))
setobject!(vis[:low], mesh_low, MeshPhongMaterial(color=Colors.RGBA(1,1,1,0.5)))
setobject!(vis[:high], mesh_high, MeshPhongMaterial(color=Colors.RGBA(0,0,0,1)))

save(joinpath(results_dir, "bluesoap_low.obj"), mesh_low)
save(joinpath(results_dir, "bluesoap_high.obj"), mesh_high)

################################################################################
# halfsoap
################################################################################
halfsoap_nerf = OSFLoader.get_nerf_object(filename="halfsoap")
mesh_low, mesh_high = nerf_mesh(halfsoap_nerf, sampling_density=100,
        normalizer=DensityFieldNormalizer(nerf=:halfsoap))
setobject!(vis[:low], mesh_low, MeshPhongMaterial(color=Colors.RGBA(1,1,1,0.5)))
setobject!(vis[:high], mesh_high, MeshPhongMaterial(color=Colors.RGBA(0,0,0,1)))

save(joinpath(results_dir, "halfsoap_low.obj"), mesh_low)
save(joinpath(results_dir, "halfsoap_high.obj"), mesh_high)
