
function build_polytope!(poly::Polyhedron, vis::Visualizer; name::Symbol=:robot, color=RGBA(1,1,1,1.0))
    mesh = Polyhedra.Mesh(poly)
    build_mesh!(mesh, vis; name=name, color=color)
end

function build_mesh!(mesh, vis::Visualizer; name::Symbol=:robot, color=RGBA(1,1,1,1.0))
    setobject!(vis[name], mesh, MeshPhongMaterial(color=color))
    return nothing
end

function set_robot!(x, q, vis::Visualizer; name::Symbol=:robot)
    settransform!(vis[name], MeshCat.compose(
        MeshCat.Translation(x),
        MeshCat.LinearMap(q),)
        )
    return nothing
end
