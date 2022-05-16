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
density = nerf_density(bluesoap_nerf, sampling_density=100,
        normalizer=DensityFieldNormalizer(nerf=:bluesoap))
for i = 0:100
    mesh = GeometryBasics.Mesh(D, NaiveSurfaceNets(iso=normalizer.density_low))
end

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



################################################################################
# Level sets
################################################################################
vis = Visualizer()
open(vis)

nerf = :bluesoap
nerf_object = OSFLoader.get_nerf_object(filename=String(nerf))
normalizer=DensityFieldNormalizer(nerf=nerf)
density = nerf_density(nerf_object, sampling_density=100, normalizer=normalizer)
levels = 0:1:100
for i in levels
    mesh = GeometryBasics.Mesh(density, NaiveSurfaceNets(iso=i))
    MeshCat.setobject!(vis[Symbol(i)], mesh)
end

anim = MeshCat.Animation(10)
for i in levels
    atframe(anim, i) do
        for j in levels
            if i == j
                setvisible!(vis[Symbol(j)], true)
            else
                setvisible!(vis[Symbol(j)], false)
            end
            settransform!(vis[Symbol(j)], MeshCat.Translation(0,0,0.6))
            set_camera!(vis, cam_pos=[2cos(pi+4π*i/length(levels)), 2sin(pi+4π*i/length(levels)), 0.5])
        end
    end
end
set_floor!(vis)
set_light!(vis)
set_background!(vis)
setanimation!(vis, anim)


convert_frames_to_video_and_gif("bunny_melting")
