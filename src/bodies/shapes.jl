"""
    Shape{T} 

    Abstract type; Subtypes contain geometric and visual information for a Body.
"""
abstract type Shape{T} end

"""
    EmptyShape{T} <: Shape{T}

    Contains no geometric or visual information
"""
struct EmptyShape{T} <: Shape{T}
    EmptyShape() = new{Float64}()
end

#TODO: change to MeshShape

"""
    Mesh{T} <: Shape{T}

    Contains geometric and visual information based on .obj file
"""
mutable struct Mesh{T} <: Shape{T}
    position_offset::SVector{3,T}
    axis_offset::Quaternion{T}

    path::String
    scale::SVector{3,T}
    color::RGBA

    function Mesh(path::String;
            position_offset::AbstractVector=szeros(3), 
            axis_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((position_offset, axis_offset))...)
        new{T}(position_offset, axis_offset, path, scale, color)
    end

    function Mesh(path::String, m::Real, J::AbstractMatrix;
            position_offset::AbstractVector=szeros(3), 
            axis_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            name::Symbol=Symbol("body_" * randstring(4)), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((m, J, position_offset, axis_offset))...)
        return Body(m, J; name=name, shape=new{T}(position_offset, axis_offset, path, scale, color))
    end
end

"""
    Box{T} <: Shape{T}

    Cuboid geometry 

    position_offset: geometry origin offset from center of mass
    axis_offset: orientation offset from body frame
    xyz: dimensions (meters)
    scale: scaling
    color: RGBA
"""
mutable struct Box{T} <: Shape{T}
    position_offset::SVector{3,T}
    axis_offset::Quaternion{T}

    xyz::SVector{3,T}
    scale::SVector{3,T}
    color::RGBA

    function Box(x::Real, y::Real, z::Real;
            position_offset::AbstractVector=szeros(3), 
            axis_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((x, y, z, position_offset, axis_offset))...)
        new{T}(position_offset, axis_offset, [x; y; z], scale, color)
    end

    function Box(x::Real, y::Real, z::Real, m::Real;
            position_offset::AbstractVector=szeros(3), 
            axis_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            name::Symbol=Symbol("body_" * randstring(4)), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((x, y, z, m, position_offset, axis_offset))...)
        J = 1 / 12 * m * diagm([y^2 + z^2; x^2 + z^2; x^2 + y^2])
        return Body(m, J; name=name, shape=new{T}(position_offset, axis_offset, [x;y;z], scale, color))
    end
end

"""
    Cylinder{T} <: Shape{T}

    cylinder geometry 
    
    position_offset: geometry origin offset from center of mass
    axis_offset: orientation offset from body frame
    rh: radius and height dimensions (meters)
    scale: scaling
    color: RGBA
"""
mutable struct Cylinder{T} <: Shape{T}
    position_offset::SVector{3,T}
    axis_offset::Quaternion{T}

    rh::SVector{2,T}
    scale::SVector{3,T}
    color::RGBA

    # Cylinder points in the z direction
    function Cylinder(r::Real, h::Real;
            position_offset::AbstractVector=szeros(3), 
            axis_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((r, h, position_offset, axis_offset))...)
        new{T}(position_offset, axis_offset, [r;h], scale, color)
    end

    function Cylinder(r::Real, h::Real, m::Real;
            position_offset::AbstractVector=szeros(3), 
            axis_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            name::Symbol=Symbol("body_" * randstring(4)), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((r, h, m, position_offset, axis_offset))...)
        J = 1 / 2 * m * diagm([r^2 + 1 / 6 * h^2; r^2 + 1 / 6 * h^2; r^2])
        return Body(m, J; name=name, shape=new{T}(position_offset, axis_offset, [r;h], scale, color))
    end
end

"""
    Capsule{T} <: Shape{T}

    capsule geometry 
    
    position_offset: geometry origin offset from center of mass
    axis_offset: orientation offset from body frame
    rh: radius and height dimensions (meters)
    scale: scaling
    color: RGBA
"""
mutable struct Capsule{T} <: Shape{T}
    position_offset::SVector{3,T}
    axis_offset::Quaternion{T}

    rh::SVector{2,T}
    scale::SVector{3,T}
    color::RGBA

    # Capsule points in the z direction
    function Capsule(r::Real, h::Real;
            position_offset::AbstractVector=szeros(3), 
            axis_offset::Quaternion= one(Quaternion),
            scale::AbstractVector=sones(3), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((r, h, position_offset, axis_offset))...)
        new{T}(position_offset, axis_offset, [r; h], scale, color)
    end

    function Capsule(r::Real, h::Real, m::Real;
            position_offset::AbstractVector=szeros(3), 
            axis_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            name::Symbol=Symbol("body_" * randstring(4)), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((r, h, m, position_offset, axis_offset))...)

        mass_cylinder = π * h * r^2.0
        mass_hemisphere = π * 2.0 / 3.0 * r^3.0 
        mass_total = mass_cylinder + 2 * mass_hemisphere
        Ixx_cylinder = mass_cylinder * (1.0 / 12.0 * h^2.0 + 1.0 / 4.0 * r^2.0)
        Ixx_hemisphere = 83.0 / 320 * mass_hemisphere * r^2
        d = (3.0 / 8.0 * r + 0.5 * h)
        Ixx = Ixx_cylinder + 2.0 * (Ixx_hemisphere + mass_hemisphere * d^2.0)
        Izz = 0.5 * mass_cylinder * r^2.0 + mass_hemisphere * (2.0 * r^2.0) / 5.0 

        J = m * diagm([Ixx; Ixx; Izz])

        return Body(m, J; name=name, shape=new{T}(position_offset, axis_offset, [r; h], scale, color))
    end
end

"""
    Shapes{T} <: Shape{T}

    composite geometry
    
    shape: list of Shape objects
    position_offset: geometry origin offset from center of mass
    axis_offset: orientation offset from body frame
    xyz: dimensions (meters)
    scale: scaling
    color: RGBA
"""
mutable struct Shapes{T} <: Shape{T}
    shape::Vector 
    position_offset::SVector{3,T}
    axis_offset::Quaternion{T}
    scale::SVector{3,T}
    color::RGBA

    function Shapes(shapes::Vector{Shape{T}}; 
        position_offset::AbstractVector=szeros(3), 
        axis_offset::Quaternion=one(Quaternion),
        scale::AbstractVector=sones(3), 
        name::Symbol=Symbol("body_" * randstring(4)), 
        color=RGBA(0.75, 0.75, 0.75)) where T
        new{T}(shapes, position_offset, axis_offset, scale, color)
    end

    function Shapes(shapes::Vector, m::T, J; 
        position_offset::AbstractVector=szeros(3), 
        axis_offset::Quaternion=one(Quaternion),
        scale::AbstractVector=sones(3), 
        name::Symbol=Symbol("body_" * randstring(4)), 
        color=RGBA(0.75, 0.75, 0.75)) where T
        Body(m, J; name=name, shape=new{T}(shapes, position_offset, axis_offset, scale, color))
    end
end

"""
    Sphere{T} <: Shape{T}

    sphere geometry 
    
    position_offset: geometry origin offset from center of mass
    axis_offset: orientation offset from body frame
    r: radius (meters)
    scale: scaling
    color: RGBA
"""
mutable struct Sphere{T} <: Shape{T}
    position_offset::SVector{3,T}
    axis_offset::Quaternion{T}
    r::T
    scale::SVector{3,T}
    color::RGBA

    function Sphere(r::Real;
            position_offset::AbstractVector=szeros(3), 
            axis_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((r, position_offset, axis_offset))...)
        new{T}(position_offset, axis_offset, r, scale, color)
    end

    function Sphere(r::Real, m::Real;
            position_offset::AbstractVector=szeros(3), 
            axis_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            name::Symbol=Symbol("body_" * randstring(4)), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((r, m, position_offset, axis_offset))...)
        J = 2 / 5 * m * diagm([r^2 for i = 1:3])
        return Body(m, J; name=name, shape=new{T}(position_offset, axis_offset, r, scale, color))
    end
end

"""
    Pyramid{T} <: Shape{T}

    pyramid geometry 
    
    position_offset: geometry origin offset from center of mass
    axis_offset: orientation offset from body frame
    wh: width and height dimensions (meters)
    scale: scaling
    color: RGBA
"""
mutable struct Pyramid{T} <: Shape{T}
    position_offset::SVector{3,T}
    axis_offset::Quaternion{T}
    wh::SVector{2,T}
    scale::SVector{3,T}
    color::RGBA

    # Pyramid points in the z direction, Center of mass at 1/4 h
    function Pyramid(w::Real, h::Real;
            position_offset::AbstractVector=szeros(3), 
            axis_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((w, h, position_offset, axis_offset))...)
        new{T}(position_offset, axis_offset, [w;h], scale, color)
    end

    function Pyramid(w::Real, h::Real, m::Real;
            position_offset::AbstractVector=szeros(3), 
            axis_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            name::Symbol=Symbol("body_" * randstring(4)), color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((w, h, m, position_offset, axis_offset))...)
        J = 1/80 * m * diagm([4*w^2+3*h^2;4*w^2+3*h^2;8*w^2])
        return Body(m, J; name=name, shape=new{T}(position_offset, axis_offset, [w;h], scale, color))
    end
end

function convert_shape(box::Box)
    x, y, z = Tuple(box.xyz)
    return GeometryBasics.HyperRectangle(Vec(-x / 2.0, -y / 2.0, -z / 2.0), Vec(x, y, z))
end

function convert_shape(cylinder::Cylinder)
    r, h = Tuple(cylinder.rh)
    return GeometryBasics.Cylinder(Point(0.0, 0.0, -h / 2.0),Point(0.0, 0.0, h / 2.0), r)
end

function convert_shape(sphere::Sphere)
    r = sphere.r
    return GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), r)
end

function convert_shape(pyramid::Pyramid)
    w, h = Tuple(pyramid.wh)
    return GeometryBasics.Pyramid(Point(0.0, 0.0, -h / 4.0), h, w)
end

function convert_shape(mesh::Mesh)
    return MeshFileObject(mesh.path)
end

function convert_shape(::EmptyShape)
    return nothing
end

function convert_shape(capsule::Capsule)
    r, h = Tuple(capsule.rh)
    p1 = Point(0.0, 0.0, -h / 2.0)
    p2 = Point(0.0, 0.0, h / 2.0)
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

