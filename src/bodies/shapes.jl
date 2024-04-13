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

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, shape::EmptyShape)
    summary(io, shape)
    println(io, "")
end

#TODO: change to MeshShape

"""
    Mesh{T} <: Shape{T}

    Contains geometric and visual information based on .obj file
"""
mutable struct Mesh{T} <: Shape{T}
    position_offset::SVector{3,T}
    orientation_offset::Quaternion{T}
    path::String
    scale::SVector{3,T}
    color::RGBA

    function Mesh(path::String;
            position_offset::AbstractVector=szeros(3), 
            orientation_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((position_offset, orientation_offset))...)
        new{T}(position_offset, orientation_offset, path, scale, color)
    end

    function Mesh(path::String, m::Real, J::AbstractMatrix;
            position_offset::AbstractVector=szeros(3), 
            orientation_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            name::Symbol=Symbol("body_" * randstring(4)), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((m, J, position_offset, orientation_offset))...)
        return Body(m, J; name=name, shape=new{T}(position_offset, orientation_offset, path, scale, color))
    end
end

"""
    Box{T} <: Shape{T}

    Cuboid geometry 

    position_offset: geometry origin offset from center of mass
    orientation_offset: orientation offset from body frame
    xyz: dimensions (meters)
    scale: scaling
    color: RGBA
"""
mutable struct Box{T} <: Shape{T}
    primitive::AbstractPrimitive
    position_offset::SVector{3,T}
    orientation_offset::Quaternion{T}
    xyz::SVector{3,T}
    scale::SVector{3,T}
    color::RGBA

    function Box(x::Real, y::Real, z::Real;
            position_offset::AbstractVector=szeros(3), 
            orientation_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((x, y, z, position_offset, orientation_offset))...)
        new{T}(box_primitive(T(x),T(y),T(z)),position_offset, orientation_offset, [x; y; z], scale, color)
    end

    function Box(x::Real, y::Real, z::Real, m::Real;
            position_offset::AbstractVector=szeros(3), 
            orientation_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            name::Symbol=Symbol("body_" * randstring(4)), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((x, y, z, m, position_offset, orientation_offset))...)
        J = 1 / 12 * m * diagm([y^2 + z^2; x^2 + z^2; x^2 + y^2])
        return Body(m, J; name=name, shape=new{T}(box_primitive(T(x),T(y),T(z)),position_offset, orientation_offset, [x;y;z], scale, color))
    end
end

"""
    Cylinder{T} <: Shape{T}

    cylinder geometry 
    
    position_offset: geometry origin offset from center of mass
    orientation_offset: orientation offset from body frame
    rh: radius and height dimensions (meters)
    scale: scaling
    color: RGBA
"""
mutable struct Cylinder{T} <: Shape{T}
    primitive::AbstractPrimitive
    position_offset::SVector{3,T}
    orientation_offset::Quaternion{T}
    rh::SVector{2,T}
    scale::SVector{3,T}
    color::RGBA

    # Cylinder points in the z direction
    function Cylinder(r::Real, h::Real;
            position_offset::AbstractVector=szeros(3), 
            orientation_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((r, h, position_offset, orientation_offset))...)
        new{T}(cylinder_primitive(T(r),T(h)),position_offset, orientation_offset, [r;h], scale, color)
    end

    function Cylinder(r::Real, h::Real, m::Real;
            position_offset::AbstractVector=szeros(3), 
            orientation_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            name::Symbol=Symbol("body_" * randstring(4)), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((r, h, m, position_offset, orientation_offset))...)
        J = 1 / 2 * m * diagm([r^2 + 1 / 6 * h^2; r^2 + 1 / 6 * h^2; r^2])
        return Body(m, J; name=name, shape=new{T}(cylinder_primitive(T(r),T(h)),position_offset, orientation_offset, [r;h], scale, color))
    end
end

"""
    Capsule geometry created as a CombinedShapes
    
    position_offset: geometry origin offset from center of mass
    orientation_offset: orientation offset from body frame
    rh: radius and height dimensions (meters)
    scale: scaling
    color: RGBA
"""
function Capsule(r::Real, h::Real;
        position_offset::AbstractVector=szeros(3), 
        orientation_offset::Quaternion=one(Quaternion),
        scale::AbstractVector=sones(3), 
        color=RGBA(0.75, 0.75, 0.75))
    T = promote_type(quateltype.((r, h, position_offset, orientation_offset))...)

    cylinder = Cylinder(r, h)
    cap1 = Sphere(r; position_offset=[0;0;h/2])
    cap2 = Sphere(r; position_offset=[0;0;-h/2])
    CombinedShapes([cylinder;cap1;cap2]; position_offset, orientation_offset, scale, color)
end

function Capsule(r::Real, h::Real, m::Real;
        position_offset::AbstractVector=szeros(3), 
        orientation_offset::Quaternion=one(Quaternion),
        scale::AbstractVector=sones(3), 
        name::Symbol=Symbol("body_" * randstring(4)), 
        color=RGBA(0.75, 0.75, 0.75))
    T = promote_type(quateltype.((r, h, m, position_offset, orientation_offset))...)

    volume_cylinder = π * h * r^2.0
    volume_hemisphere = π * 4.0 / 3.0 * r^3.0 / 2.0
    volume_total = volume_cylinder + 2 * volume_hemisphere
    mass_cylinder = m * volume_cylinder / volume_total
    mass_hemisphere = m * volume_hemisphere / volume_total
    Ixx_cylinder = mass_cylinder * (1.0 / 12.0 * h^2.0 + 1.0 / 4.0 * r^2.0)
    Izz_cylinder = mass_cylinder * 1.0 / 2.0 * r^2.0
    Ixx_hemisphere = 83.0 / 320 * mass_hemisphere * r^2
    Izz_hemisphere = mass_hemisphere * 2.0 / 5.0 * r^2 / 2.0
    d = (3.0 / 8.0 * r + 0.5 * h)
    Ixx = Ixx_cylinder + 2.0 * (Ixx_hemisphere + mass_hemisphere * d^2.0)
    Izz = Izz_cylinder + Izz_hemisphere * 2.0

    J = diagm([Ixx; Ixx; Izz])

    return Body(m, J; name=name, shape=Capsule(r, h; position_offset, orientation_offset, scale, color))
end

"""
CombinedShapes{T} <: Shape{T}

    composite geometry
    
    position_offset: geometry origin offset from center of mass
    orientation_offset: orientation offset from body frame
    shape: list of Shape objects
    xyz: dimensions (meters)
"""
mutable struct CombinedShapes{T} <: Shape{T}
    position_offset::SVector{3,T}
    orientation_offset::Quaternion{T}
    shapes::Vector{<:Shape{T}}
    scale::SVector{3,T}
    color::RGBA

    function CombinedShapes(shapes::Vector{<:Shape{T}}; 
        position_offset::AbstractVector=szeros(3), 
        orientation_offset::Quaternion=one(Quaternion),
        scale::AbstractVector=sones(3), 
        color=RGBA(0.75, 0.75, 0.75)) where T

        new{T}(position_offset, orientation_offset, shapes, scale, color)
    end

    function CombinedShapes(shapes::Vector{<:Shape{T}}, m::T, J; 
        position_offset::AbstractVector=szeros(3), 
        orientation_offset::Quaternion=one(Quaternion),
        name::Symbol=Symbol("body_" * randstring(4)),
        scale::AbstractVector=sones(3), 
        color=RGBA(0.75, 0.75, 0.75)) where T

        Body(m, J; name=name, shape=new{T}(position_offset, orientation_offset, shapes, scale, color))
    end
end

"""
    Sphere{T} <: Shape{T}

    sphere geometry 
    
    position_offset: geometry origin offset from center of mass
    orientation_offset: orientation offset from body frame
    r: radius (meters)
    scale: scaling
    color: RGBA
"""
mutable struct Sphere{T} <: Shape{T}
    primitive::AbstractPrimitive
    position_offset::SVector{3,T}
    orientation_offset::Quaternion{T}
    r::T
    scale::SVector{3,T}
    color::RGBA

    function Sphere(r::Real;
            position_offset::AbstractVector=szeros(3), 
            orientation_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((r, position_offset, orientation_offset))...)
        new{T}(sphere_primitive(T(r)),position_offset, orientation_offset, r, scale, color)
    end

    function Sphere(r::Real, m::Real;
            position_offset::AbstractVector=szeros(3), 
            orientation_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            name::Symbol=Symbol("body_" * randstring(4)), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((r, m, position_offset, orientation_offset))...)
        J = 2 / 5 * m * diagm([r^2 for i = 1:3])
        return Body(m, J; name=name, shape=new{T}(sphere_primitive(T(r)),position_offset, orientation_offset, r, scale, color))
    end
end

"""
    Pyramid{T} <: Shape{T}

    pyramid geometry 
    
    position_offset: geometry origin offset from center of mass
    orientation_offset: orientation offset from body frame
    wh: width and height dimensions (meters)
    scale: scaling
    color: RGBA
"""
mutable struct Pyramid{T} <: Shape{T}
    primitive::AbstractPrimitive
    position_offset::SVector{3,T}
    orientation_offset::Quaternion{T}
    wh::SVector{2,T}
    scale::SVector{3,T}
    color::RGBA

    # Pyramid points in the z direction, Center of mass at 1/4 h
    function Pyramid(w::Real, h::Real;
            position_offset::AbstractVector=szeros(3), 
            orientation_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((w, h, position_offset, orientation_offset))...)
        new{T}(pyramid_primitive(T(w),T(h)),position_offset, orientation_offset, [w;h], scale, color)
    end

    function Pyramid(w::Real, h::Real, m::Real;
            position_offset::AbstractVector=szeros(3), 
            orientation_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            name::Symbol=Symbol("body_" * randstring(4)), color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((w, h, m, position_offset, orientation_offset))...)
        J = 1/80 * m * diagm([4*w^2+3*h^2;4*w^2+3*h^2;8*w^2])
        return Body(m, J; name=name, shape=new{T}(pyramid_primitive(T(w),T(h)),position_offset, orientation_offset, [w;h], scale, color))
    end
end

"""
    FrameShape{T} <: Shape{T}

    coordinate frame geometry 
    
    position_offset: geometry origin offset from center of mass
    orientation_offset: orientation offset from body frame
    scale: scaling
    color: not used
"""
mutable struct FrameShape{T} <: Shape{T}
    position_offset::SVector{3,T}
    orientation_offset::Quaternion{T}
    scale::SVector{3,T}
    color::RGBA

    function FrameShape(;
            position_offset::AbstractVector=szeros(3), 
            orientation_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((position_offset, orientation_offset))...)
        new{T}(position_offset, orientation_offset, scale, color)
    end

    function FrameShape(m::Real;
            position_offset::AbstractVector=szeros(3), 
            orientation_offset::Quaternion=one(Quaternion),
            scale::AbstractVector=sones(3), 
            name::Symbol=Symbol("body_" * randstring(4)), color=RGBA(0.75, 0.75, 0.75))
        T = promote_type(quateltype.((m, position_offset, orientation_offset))...)
        J = m * sI(3)
        return Body(m, J; name=name, shape=new{T}(position_offset, orientation_offset, scale, color))
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, shape::Shape)
    summary(io, shape)
    println(io,"")
    println(io," position_offset: "*string(shape.position_offset))
    println(io," orientation_offset:     "*string(shape.orientation_offset))
    println(io," scale:           "*string(shape.scale))
    println(io," color:           "*string(shape.color))
end

# TODO position and orientation offset not accounted for
function box_primitive(x,y,z)
    A = SA[
        1 0 0
        0 1 0 
        0 0 1
        -1 0 0
        0 -1 0
        0 0 -1.0
    ]
    b = SA[x;y;z;x;y;z]/2

    return DifferentiableCollisions.Polytope(A,b)
end

function cylinder_primitive(r,h)
    return DifferentiableCollisions.Cylinder(r, h)
end

function sphere_primitive(r)
    return DifferentiableCollisions.Sphere(r)
end

function pyramid_primitive(w,h)
    A = SA[
        1/(3/8*w)  0          1/(3/4*h)
        0          1/(3/8*w)  1/(3/4*h)
       -1/(3/8*w)  0          1/(3/4*h)
        0         -1/(3/8*w)  1/(3/4*h)
        0          0         -1/(1/4*h)
    ]
    b = SA[1;1;1;1;1.0]

    return DifferentiableCollisions.Polytope(A,b)
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

function convert_shape(frame::FrameShape)
    return MeshCat.Triad()
end

function convert_shape(mesh::Mesh)
    return MeshFileObject(mesh.path)
end

function convert_shape(::EmptyShape)
    return nothing
end

function convert_shape(combinedshapes::CombinedShapes)
    geom = []
    for shape in combinedshapes.shapes
        push!(geom, convert_shape(shape))
    end
    return geom
end

# color
set_color!(shape::EmptyShape, color) = nothing

function set_color!(shape::Shape, color)
    shape.color = color
end

function set_color!(shapes::CombinedShapes, color)
    for i in eachindex(shapes)
        shapes.shape[i].color = color
    end
end

