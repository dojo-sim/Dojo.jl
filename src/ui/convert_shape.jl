function convert_shape(box::Box)
    x,y,z = Tuple(box.xyz)
    return GeometryBasics.HyperRectangle(Vec(-x/2,-y/2,-z/2),Vec(x,y,z))
end

function convert_shape(cylinder::Cylinder)
    r,h = Tuple(cylinder.rh)
    return GeometryBasics.Cylinder(Point(0.0,0.0,-h/2),Point(0.0,0.0,h/2), r)
end

function convert_shape(sphere::Sphere)
    r = sphere.r
    return GeometryBasics.Sphere(Point(0.0,0.0,0.0), r)
end

function convert_shape(pyramid::Pyramid)
    w, h = Tuple(pyramid.wh)
    return GeometryBasics.Pyramid(Point(0.0,0.0,-h/4), h, w)
end

function convert_shape(mesh::Mesh)
    return MeshFileObject(joinpath(@__DIR__, "..", "..", mesh.path))
end

function convert_shape(::EmptyShape)
    return nothing
end

function convert_shape(capsule::Capsule)
    r, h = Tuple(capsule.rh)
    p1 = Point(0.0, 0.0, -h / 2)
    p2 = Point(0.0, 0.0, h / 2)
    cyl = GeometryBasics.Cylinder(p1, p2, r)
    cap1 = GeometryBasics.Sphere(p1, r)
    cap2 = GeometryBasics.Sphere(p2, r)
    return [cyl, cap1, cap2]
end

function convert_shape(shapes::Shapes)
    geom = []
    for s in shapes.shape
        push!(geom, convert_shape(s))
    end
    return geom
end
