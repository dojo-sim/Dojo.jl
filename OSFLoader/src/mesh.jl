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

function nerf_mesh(nerf_object;
    normalizer::DensityFieldNormalizer=DensityFieldNormalizer(),
    sampling_density::Int=20)

    D = nerf_density(nerf_object, normalizer=normalizer, sampling_density=sampling_density)
    mesh_low = GeometryBasics.Mesh(D, NaiveSurfaceNets(iso=normalizer.density_low))
    mesh_high = GeometryBasics.Mesh(D, NaiveSurfaceNets(iso=normalizer.density_high))

    return mesh_low, mesh_high
end
