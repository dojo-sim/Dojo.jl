using OSFLoader
using GeometryBasics
using LinearAlgebra
using JLD2
using MeshIO
using Meshing
using MeshCat
using Quaternions
using MeshIO
using Polyhedra
using GLPK

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
for i = 0:20:100
    mesh = GeometryBasics.Mesh(density, NaiveSurfaceNets(iso=i))
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

nerf = :bunny
nerf_object = OSFLoader.get_nerf_object(filename=String(nerf))
normalizer = DensityFieldNormalizer(nerf=nerf)
density = nerf_density(nerf_object, sampling_density=100, normalizer=normalizer)
levels = [0:1:10; 20:10:300]
for i in levels
    mesh = GeometryBasics.Mesh(density, MarchingCubes(iso=i))
    MeshCat.setobject!(vis[Symbol(i)], mesh)
    MeshIO.save(joinpath(@__DIR__, "bunny_level_$i.obj"), mesh)
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
            # Dojo.set_camera!(vis, cam_pos=[2cos(pi+4π*i/length(levels)), 2sin(pi+4π*i/length(levels)), 0.5])
        end
    end
end
Dojo.set_floor!(vis)
Dojo.set_light!(vis)
Dojo.set_background!(vis)
setanimation!(vis, anim)


################################################################################
# Intersection over union
################################################################################
vis = Visualizer()
render(vis)

function level_set_volume(nerf::Symbol, levels; sampling_density::Int=50)
    volumes = zeros(length(levels))

    nerf_object = OSFLoader.get_nerf_object(filename=String(nerf))
    normalizer = DensityFieldNormalizer(nerf=nerf)
    density = nerf_density(nerf_object, sampling_density=sampling_density,
        normalizer=normalizer)

    for (i,level) in enumerate(levels)
        volumes[i] = sum(density .>= level) / sampling_density^3
    end
    return volumes
end

function convex_hull_volume(nerf::Symbol, levels; sampling_density::Int=50)
    volumes = zeros(length(levels))

    for (l,level) in enumerate(levels)
        hull = convex_hull(nerf, level, sampling_density=sampling_density)
        for i = 1:sampling_density
            for j = 1:sampling_density
                for k = 1:sampling_density
                    point = (([i,j,k] .- 1) ./ (sampling_density-1))
                    (point ∈ hull) && (volumes[l] += 1 / sampling_density^3)
                end
            end
        end
    end
    return volumes
end

function convex_hull(nerf::Symbol, level; sampling_density::Int=50)
    nerf_object = OSFLoader.get_nerf_object(filename=String(nerf))
    normalizer = DensityFieldNormalizer(nerf=nerf)
    density = nerf_density(nerf_object, sampling_density=sampling_density,
        normalizer=normalizer)
    points = []
    for i = 1:sampling_density
        for j = 1:sampling_density
            for k = 1:sampling_density
                (density[k,j,i] > level) && push!(points, ([k,j,i] .-1) ./ (sampling_density-1))
            end
        end
    end
    v = convexhull(points...)
    @warn "computing convex hull..."
    v = removevredundancy(v, GLPK.Optimizer)
    p = polyhedron(v)
    hrep(p)
    return p
end

function level_set_iou(nerf::Symbol, levels, hull; sampling_density::Int=50)
    intersection = zeros(length(levels))
    union = zeros(length(levels))

    nerf_object = OSFLoader.get_nerf_object(filename=String(nerf))
    normalizer = DensityFieldNormalizer(nerf=nerf)
    density = nerf_density(nerf_object, sampling_density=sampling_density,
        normalizer=normalizer)

    for (l,level) in enumerate(levels)
        for i = 1:sampling_density
            for j = 1:sampling_density
                for k = 1:sampling_density
                    point = (([i,j,k] .- 1) ./ (sampling_density-1))
                    if point ∈ hull && density[k,j,i] >= level
                        intersection[l] += 1 / sampling_density^3
                    end
                    if point ∈ hull || density[k,j,i] >= level
                        union[l] += 1 / sampling_density^3
                    end
                end
            end
        end
    end
    return intersection ./ union, intersection, union
end

function convex_hull_iou(nerf::Symbol, levels, hull; sampling_density::Int=50)
    intersection = zeros(length(levels))
    union = zeros(length(levels))

    for (l,level) in enumerate(levels)
        level_hull = convex_hull(nerf, level, sampling_density=sampling_density)
        for i = 1:sampling_density
            for j = 1:sampling_density
                for k = 1:sampling_density
                    point = (([i,j,k] .- 1) ./ (sampling_density-1))
                    if point ∈ hull && point ∈ level_hull
                        intersection[l] += 1 / sampling_density^3
                    end
                    if point ∈ hull || point ∈ level_hull
                        union[l] += 1 / sampling_density^3
                    end
                end
            end
        end
    end
    return intersection ./ union, intersection, union
end


nerf = :bluesoap
nerf_object = OSFLoader.get_nerf_object(filename=String(nerf))
normalizer = DensityFieldNormalizer(nerf=nerf)
density = nerf_density(nerf_object, sampling_density=50, normalizer=normalizer)
levels = [0:1:10; 20:10:200]
for i in reverse(levels)
    mesh = GeometryBasics.Mesh(density, MarchingCubes(iso=i))
    MeshCat.setobject!(vis[Symbol(i)], mesh)
    # MeshIO.save(joinpath(@__DIR__, "bunny_level_$i.obj"), mesh)
end



p = convex_hull(nerf, 40.0)
mesh = Polyhedra.Mesh(p)
MeshCat.setobject!(vis[:hull], mesh)

levels = [0:1:10; 20:10:200]
normalized_volumes = normalized_volume(nerf, levels, sampling_density=50)
scatter(levels, normalized_volumes)

levels = 50:10:200
cvx_volumes = convex_hull_volume(nerf, levels)
scatter(levels, cvx_volumes)



# levels = [10:10:10; 20:10:200]
# open(vis)
# for i in reverse(levels)
#     mesh = Polyhedra.Mesh(convex_hull(nerf, i))
#     MeshCat.setobject!(vis[Symbol(i)], mesh)
#     # MeshIO.save(joinpath(@__DIR__, "bunny_level_$i.obj"), mesh)
# end


levels = [1:2:10; 20:10:200]
hull70 = convex_hull(nerf, 70)
hull65 = convex_hull(nerf, 65)
lvl_iou, lvl_i, lvl_u = level_set_iou(nerf, levels, hull70)
scatter(levels, lvl_iou)
scatter!(levels, lvl_i)
scatter!(levels, lvl_u)
cvx_iou, cvx_i, cvx_u = convex_hull_iou(nerf, levels, hull65)
scatter(levels, cvx_iou)
scatter!(levels, cvx_i)
scatter!(levels, cvx_u)
