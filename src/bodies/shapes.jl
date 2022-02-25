abstract type Shape{T} end

struct EmptyShape{T} <: Shape{T}
    EmptyShape() = new{Float64}()
end

mutable struct Mesh{T} <: Shape{T}
    xoffset::SVector{3,T}
    axis_offset::UnitQuaternion{T}

    path::String
    scale::SVector{3,T}
    color::RGBA

    function Mesh(path::String;
            xoffset::AbstractVector=szeros(3), axis_offset::UnitQuaternion=one(UnitQuaternion),
            scale::AbstractVector=sones(3), color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(eltype.((xoffset, axis_offset))...)
        new{T}(xoffset, axis_offset, path, scale, color)
    end

    function Mesh(path::String, m::Real, J::AbstractMatrix;
            xoffset::AbstractVector=szeros(3), axis_offset::UnitQuaternion=one(UnitQuaternion),
            scale::AbstractVector=sones(3), name::Symbol=Symbol("body_" * randstring(4)), color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(eltype.((m, J, xoffset, axis_offset))...)
        return Body(m, J; name=name, shape=new{T}(xoffset, axis_offset, path, scale, color))
    end
end

mutable struct Box{T} <: Shape{T}
    xoffset::SVector{3,T}
    axis_offset::UnitQuaternion{T}

    xyz::SVector{3,T}
    scale::SVector{3,T}
    color::RGBA

    function Box(x::Real, y::Real, z::Real;
            xoffset::AbstractVector=szeros(3), axis_offset::UnitQuaternion=one(UnitQuaternion),
            scale::AbstractVector=sones(3), color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(eltype.((x, y, z, xoffset, axis_offset))...)
        new{T}(xoffset, axis_offset, [x;y;z], scale, color)
    end

    function Box(x::Real, y::Real, z::Real, m::Real;
            xoffset::AbstractVector=szeros(3), axis_offset::UnitQuaternion=one(UnitQuaternion),
            scale::AbstractVector=sones(3), name::Symbol=Symbol("body_" * randstring(4)), color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(eltype.((x, y, z, m, xoffset, axis_offset))...)
        J = 1 / 12 * m * diagm([y^2 + z^2; x^2 + z^2; x^2 + y^2])
        return Body(m, J; name=name, shape=new{T}(xoffset, axis_offset, [x;y;z], scale, color))
    end
end

mutable struct Cylinder{T} <: Shape{T}
    xoffset::SVector{3,T}
    axis_offset::UnitQuaternion{T}

    rh::SVector{2,T}
    scale::SVector{3,T}
    color::RGBA

    # Cylinder points in the z direction
    function Cylinder(r::Real, h::Real;
            xoffset::AbstractVector=szeros(3), axis_offset::UnitQuaternion=one(UnitQuaternion),
            scale::AbstractVector=sones(3), color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(eltype.((r, h, xoffset, axis_offset))...)
        new{T}(xoffset, axis_offset, [r;h], scale, color)
    end

    function Cylinder(r::Real, h::Real, m::Real;
            xoffset::AbstractVector=szeros(3), axis_offset::UnitQuaternion=one(UnitQuaternion),
            scale::AbstractVector=sones(3), name::Symbol=Symbol("body_" * randstring(4)), color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(eltype.((r, h, m, xoffset, axis_offset))...)
        J = 1 / 2 * m * diagm([r^2 + 1 / 6 * h^2; r^2 + 1 / 6 * h^2; r^2])
        return Body(m, J; name=name, shape=new{T}(xoffset, axis_offset, [r;h], scale, color))
    end
end

mutable struct Capsule{T} <: Shape{T}
    xoffset::SVector{3,T}
    axis_offset::UnitQuaternion{T}

    rh::SVector{2,T}
    scale::SVector{3,T}
    color::RGBA

    # Capsule points in the z direction
    function Capsule(r::Real, h::Real;
            xoffset::AbstractVector=szeros(3), axis_offset::UnitQuaternion= one(UnitQuaternion),
            scale::AbstractVector=sones(3), color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(eltype.((r, h, xoffset, axis_offset))...)
        new{T}(xoffset, axis_offset, [r; h], scale, color)
    end

    function Capsule(r::Real, h::Real, m::Real;
            xoffset::AbstractVector=szeros(3), axis_offset::UnitQuaternion=one(UnitQuaternion),
            scale::AbstractVector=sones(3), name::Symbol=Symbol("body_" * randstring(4)), color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(eltype.((r, h, m, xoffset, axis_offset))...)

        mass_cylinder = π * h * r^2.0
        mass_hemisphere = π * 2.0 / 3.0 * r^3.0 
        mass_total = mass_cylinder + 2 * mass_hemisphere
        Ixx_cylinder = mass_cylinder * (1.0 / 12.0 * h^2.0 + 1.0 / 4.0 * r^2.0)
        Ixx_hemisphere = 83.0 / 320 * mass_hemisphere * r^2
        d = (3.0 / 8.0 * r + 0.5 * h)
        Ixx = Ixx_cylinder + 2.0 * (Ixx_hemisphere + mass_hemisphere * d^2.0)
        Izz = 0.5 * mass_cylinder * r^2.0 + mass_hemisphere * (2.0 * r^2.0) / 5.0 

        J = m * diagm([Ixx; Ixx; Izz])

        return Body(m, J; name=name, shape=new{T}(xoffset, axis_offset, [r; h], scale, color))
    end
end

mutable struct Shapes{T} <: Shape{T}
    shape::Vector 
    xoffset::SVector{3,T}
    axis_offset::UnitQuaternion{T}
    scale::SVector{3,T}
    color::RGBA

    function Shapes(shapes::Vector{Shape{T}}; 
        xoffset::AbstractVector=szeros(3), axis_offset::UnitQuaternion=one(UnitQuaternion),
        scale::AbstractVector=sones(3), name::Symbol=Symbol("body_" * randstring(4)), color=RGBA(0.75, 0.75, 0.75)) where T
        new{T}(shapes, xoffset, axis_offset, scale, color)
    end

    function Shapes(shapes::Vector, m::T, J; 
        xoffset::AbstractVector=szeros(3), axis_offset::UnitQuaternion=one(UnitQuaternion),
        scale::AbstractVector=sones(3), name::Symbol=Symbol("body_" * randstring(4)), color=RGBA(0.75, 0.75, 0.75)) where T
        Body(m, J; name=name, shape=new{T}(shapes, xoffset, axis_offset, scale, color))
    end
end

mutable struct Sphere{T} <: Shape{T}
    xoffset::SVector{3,T}
    axis_offset::UnitQuaternion{T}

    r::T
    scale::SVector{3,T}
    color::RGBA

    function Sphere(r::Real;
            xoffset::AbstractVector=szeros(3), axis_offset::UnitQuaternion=one(UnitQuaternion),
            scale::AbstractVector=sones(3), color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(eltype.((r, xoffset, axis_offset))...)
        new{T}(xoffset, axis_offset, r, scale, color)
    end

    function Sphere(r::Real, m::Real;
            xoffset::AbstractVector=szeros(3), axis_offset::UnitQuaternion=one(UnitQuaternion),
            scale::AbstractVector=sones(3), name::Symbol=Symbol("body_" * randstring(4)), color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(eltype.((r, m, xoffset, axis_offset))...)
        J = 2 / 5 * m * diagm([r^2 for i = 1:3])
        return Body(m, J; name=name, shape=new{T}(xoffset, axis_offset, r, scale, color))
    end
end

mutable struct Pyramid{T} <: Shape{T}
    xoffset::SVector{3,T}
    axis_offset::UnitQuaternion{T}

    wh::SVector{2,T}
    scale::SVector{3,T}
    color::RGBA

    # Pyramid points in the z direction, Center of mass at 1/4 h
    function Pyramid(w::Real, h::Real;
            xoffset::AbstractVector=szeros(3), axis_offset::UnitQuaternion=one(UnitQuaternion),
            scale::AbstractVector=sones(3), color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(eltype.((w, h, xoffset, axis_offset))...)
        new{T}(xoffset, axis_offset, [w;h], scale, color)
    end

    function Pyramid(w::Real, h::Real, m::Real;
            xoffset::AbstractVector=szeros(3), axis_offset::UnitQuaternion=one(UnitQuaternion),
            scale::AbstractVector=sones(3), name::Symbol=Symbol("body_" * randstring(4)), color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(eltype.((w, h, m, xoffset, axis_offset))...)
        J = 1/80 * m * diagm([4*w^2+3*h^2;4*w^2+3*h^2;8*w^2])
        return Body(m, J; name=name, shape=new{T}(xoffset, axis_offset, [w;h], scale, color))
    end
end

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

# color
set_color!(shape::EmptyShape, color) = nothing

function set_color!(shape::Shape, color)
    shape.color = color
end

function set_color!(shapes::Shapes, color)
    for i in eachindex(shapes)
        shapes.shape[i].color = color
    end
end

