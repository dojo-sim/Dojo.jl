using Dojo
using OSFLoader
using PyCall
using GeometryBasics
using LinearAlgebra
using JLD2
using MeshIO


function MeshIO.save(f::MeshIO.Stream{MeshIO.format"OBJ"}, mesh::AbstractMesh)
    io = MeshIO.stream(f)
    for p in decompose(Point3f0, mesh)
        println(io, "v ", p[1], " ", p[2], " ", p[3])
    end

    if hasproperty(mesh, :uv)
        for uv in mesh.uv
            println(io, "vt ", uv[1], " ", uv[2])
        end
    end

    if hasproperty(mesh, :normals)
        for n in decompose_normals(mesh)
            println(io, "vn ", n[1], " ", n[2], " ", n[3])
        end
    end

    F = eltype(faces(mesh))
    for f in decompose(F, mesh)
        println(io, "f ", join(convert.(Int, f), " "))
    end
end


vis = Visualizer()
open(vis)

################################################################################
# bunny
################################################################################
bunny_nerf = OSFLoader.get_nerf_object(filename="bunny_trans")

# NeRF density -> density tensor -> Nesh
N = 100
xrange = range(-1.0, stop=1.0, length=N)
yrange = range(-1.0, stop=1.0, length=N)
zrange = range(-1.0, stop=1.0, length=N)
D = zeros(N,N,N)
for i = 1:N
    D[:,:,i] = slice_density(bunny_nerf, xrange, yrange, zrange[i])
end

mesh = GeometryBasics.Mesh(D, NaiveSurfaceNets(iso=65))
setobject!(vis[:old], mesh,
    MeshPhongMaterial(color=RGBA{Float32}(1, 0, 0, 0.5)))


################################################################################
# bluesoap
################################################################################
bluesoap_nerf = OSFLoader.get_nerf_object(filename="bluesoap")

# NeRF density -> density tensor -> Nesh
N = 100
p = [-0.0207,  2.5801,  4.9728]
xrange = p[1] .+ range(-2.5, stop=2.5, length=N)
yrange = p[2] .+ range(-2.5, stop=2.5, length=N)
zrange = p[3] .+ range(-2.5, stop=2.5, length=N)
D = zeros(N,N,N)
for i = 1:N
    D[:,:,i] = slice_density(bluesoap_nerf, xrange, yrange, zrange[i])
end

for k in 1:111
    @show k
    mesh = GeometryBasics.Mesh(D, NaiveSurfaceNets(iso=k))
    setobject!(vis[:old], mesh,
        MeshPhongMaterial(color=RGBA{Float32}(1, 0, 0, 0.5)))
end

mesh = GeometryBasics.Mesh(D, NaiveSurfaceNets(iso=24))
tmpdir = @__DIR__
save(joinpath(tmpdir, "bluesoap.obj"), mesh)
mesh

rect = GeometryBasics.Rect(MeshCat.Vec(p...), MeshCat.Vec(5,5,5.0))

function point_density(x::Point, nerf_object)
    x = [x[1] x[2] x[3]]
    x = convert.(Float32, x)
    d = density_query(nerf_object, x)
    return d[1]
end

point_density(Point(1,2,3.0), bluesoap_nerf)

mesh = GeometryBasics.Mesh(x -> point_density(x, bluesoap_nerf), rect, MarchingCubes(iso=24), samples=(30,30,30))
setobject!(vis[:old], mesh,
    MeshPhongMaterial(color=RGBA{Float32}(1, 0, 0, 0.5)))


################################################################################
# halfsoap
################################################################################
halfsoap_nerf = OSFLoader.get_nerf_object(filename="halfsoap")

# NeRF density -> density tensor -> Nesh
N = 100
p = [-2.5519,  2.4537,  4.8040]
xrange = p[1] .+ range(-2.5, stop=2.5, length=N)
yrange = p[2] .+ range(-2.5, stop=2.5, length=N)
zrange = p[3] .+ range(-2.5, stop=2.5, length=N)
D = zeros(N,N,N)
for i = 1:N
    D[:,:,i] = slice_density(halfsoap_nerf, xrange, yrange, zrange[i])
end

for k in -100:1
    @show k
    mesh = GeometryBasics.Mesh(D, NaiveSurfaceNets(iso=k))
    setobject!(vis[:old], mesh,
        MeshPhongMaterial(color=RGBA{Float32}(1, 0, 0, 0.5)))
end


mesh = GeometryBasics.Mesh(D, NaiveSurfaceNets(iso=4))
tmpdir = @__DIR__
save(joinpath(tmpdir, "halfsoap.obj"), mesh)
mesh
