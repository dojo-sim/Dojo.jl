function convertshape(box::Box)
    x,y,z = Tuple(box.xyz)
    return GeometryBasics.HyperRectangle(Vec(-x/2,-y/2,-z/2),Vec(x,y,z))
end

function convertshape(cylinder::Cylinder)
    r,h = Tuple(cylinder.rh)
    return GeometryBasics.Cylinder(Point(0.0,0.0,-h/2),Point(0.0,0.0,h/2), r)
end

function convertshape(sphere::Sphere)
    r = sphere.r
    return GeometryBasics.Sphere(Point(0.0,0.0,0.0), r)
end

function convertshape(pyramid::Pyramid)
    w, h = Tuple(pyramid.wh)
    return GeometryBasics.Pyramid(Point(0.0,0.0,-h/4), h, w)
end

function convertshape(mesh::Mesh)
    return MeshFileObject(joinpath(@__DIR__, "..", "..", mesh.path))
end

function convertshape(::EmptyShape)
    return nothing
end

function convertshape(capsule::Capsule) 
    r, h = Tuple(capsule.rh)
    p1 = Point(0.0, 0.0, -h / 2)
    p2 = Point(0.0, 0.0, h / 2)
    cyl = GeometryBasics.Cylinder(p1, p2, r)
    cap1 = GeometryBasics.Sphere(p1, r)
    cap2 = GeometryBasics.Sphere(p2, r)
    return [cyl, cap1, cap2] 
end

function convertshape(shapes::Shapes14) 
    geom = []
    for s in shapes.shape
        push!(geom, convertshape(s))
    end
    return geom
end
