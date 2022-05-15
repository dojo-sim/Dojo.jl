using GeometryBasics
using MeshIO
import MeshIO.save
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

outer_mesh = jldopen(joinpath(example_dir(), "results/bunny_outer_mesh.jld2"))["mesh"]
inner_mesh = jldopen(joinpath(example_dir(), "results/bunny_inner_mesh.jld2"))["mesh"]

filepath_outer = joinpath(module_dir(), "environments/bunny/deps/bunny_outer_mesh.obj")
filepath_inner = joinpath(module_dir(), "environments/bunny/deps/bunny_inner_mesh.obj")
MeshIO.save(MeshIO.File{MeshIO.format"OBJ"}(filepath_outer), outer_mesh)
MeshIO.save(MeshIO.File{MeshIO.format"OBJ"}(filepath_inner), inner_mesh)
