"""
$(TYPEDEF)

A `Shape` is used for visualization of a [`Body`](@ref). All `Shape`s have an `xoffset` and a `qoffset` parameter to specify the position and orientation offset between the actual [`Body`](@ref) and the visualized `Shape`.
# Important attributes
* `xoffset`: The relative position offset (vector) between `Body` and `Shape`.  
* `qoffset`: The relative orientation offset (quaternion) between `Body` and `Shape`.
* `color`:   The color of the `Shape`.
"""
abstract type Shape{T} end

struct EmptyShape{T} <: Shape{T}
    EmptyShape() = new{Float64}()
end

"""
$(TYPEDEF)

A `Mesh` can be used to visualize arbitrary geometries.
# Important attributes
* `path`: The path to the geometry file.  
"""
mutable struct Mesh{T} <: Shape{T}
    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    path::String
    scale::SVector{3,T}
    color::RGBA

    function Mesh(path::String;
            xoffset::AbstractVector = szeros(3), qoffset::UnitQuaternion = one(UnitQuaternion),
            scale::AbstractVector = sones(3), color = RGBA(0.75, 0.75, 0.75)
        )
        T = promote_type(eltype.((xoffset, qoffset))...)

        new{T}(xoffset, qoffset, path, scale, color)
    end

    function Mesh(path::String, m::Real, J::AbstractMatrix;
            xoffset::AbstractVector = szeros(3), qoffset::UnitQuaternion = one(UnitQuaternion),
            scale::AbstractVector = sones(3), name::String="", color = RGBA(0.75, 0.75, 0.75)
        )
        T = promote_type(eltype.((m, J, xoffset, qoffset))...)

        return Body(m, J; name=name, shape=new{T}(xoffset, qoffset, path, scale, color))
    end
end

"""
$(TYPEDEF)

A `Box`.
# Important attributes
* `xyz`: The box size in x, y, and z direction (vector).  
"""
mutable struct Box{T} <: Shape{T}
    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    xyz::SVector{3,T}
    scale::SVector{3,T}
    color::RGBA


    function Box(x::Real, y::Real, z::Real;
            xoffset::AbstractVector = szeros(3), qoffset::UnitQuaternion = one(UnitQuaternion),
            scale::AbstractVector = sones(3), color = RGBA(0.75, 0.75, 0.75)
        )
        T = promote_type(eltype.((x, y, z, xoffset, qoffset))...)

        new{T}(xoffset, qoffset, [x;y;z], scale, color)
    end

    function Box(x::Real, y::Real, z::Real, m::Real;
            xoffset::AbstractVector = szeros(3), qoffset::UnitQuaternion = one(UnitQuaternion),
            scale::AbstractVector = sones(3), name::String="", color = RGBA(0.75, 0.75, 0.75)
        )
        T = promote_type(eltype.((x, y, z, m, xoffset, qoffset))...)
        J = 1 / 12 * m * diagm([y^2 + z^2;x^2 + z^2;x^2 + y^2])

        return Body(m, J; name=name, shape=new{T}(xoffset, qoffset, [x;y;z], scale, color))
    end
end

"""
$(TYPEDEF)

A `Cylinder` along the z-axis.
# Important attributes
* `rh`: The cylinder radius and height (vector).  
"""
mutable struct Cylinder{T} <: Shape{T}
    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    rh::SVector{2,T}
    scale::SVector{3,T}
    color::RGBA

    # Cylinder points in the z direction
    function Cylinder(r::Real, h::Real;
            xoffset::AbstractVector = szeros(3), qoffset::UnitQuaternion = one(UnitQuaternion),
            scale::AbstractVector = sones(3), color = RGBA(0.75, 0.75, 0.75)
        )
        T = promote_type(eltype.((r, h, xoffset, qoffset))...)

        new{T}(xoffset, qoffset, [r;h], scale, color)
    end

    function Cylinder(r::Real, h::Real, m::Real;
            xoffset::AbstractVector = szeros(3), qoffset::UnitQuaternion = one(UnitQuaternion),
            scale::AbstractVector = sones(3), name::String="", color = RGBA(0.75, 0.75, 0.75)
        )
        T = promote_type(eltype.((r, h, m, xoffset, qoffset))...)
        J = 1 / 2 * m * diagm([r^2 + 1 / 6 * h^2;r^2 + 1 / 6 * h^2;r^2])

        return Body(m, J; name=name, shape=new{T}(xoffset, qoffset, [r;h], scale, color))
    end
end

"""
$(TYPEDEF)

A `Capsule` along the z-axis.
# Important attributes
* `rh`: The cylinder radius and height (vector)
"""
mutable struct Capsule{T} <: Shape{T}
    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    rh::SVector{2,T}
    scale::SVector{3,T}
    color::RGBA

    # Capsule points in the z direction
    function Capsule(r::Real, h::Real;
            xoffset::AbstractVector = szeros(3), qoffset::UnitQuaternion = one(UnitQuaternion),
            scale::AbstractVector = sones(3), color = RGBA(0.75, 0.75, 0.75)
        )
        T = promote_type(eltype.((r, h, xoffset, qoffset))...)

        new{T}(xoffset, qoffset, [r; h], scale, color)
    end

    function Capsule(r::Real, h::Real, m::Real;
            xoffset::AbstractVector = szeros(3), qoffset::UnitQuaternion = one(UnitQuaternion),
            scale::AbstractVector = sones(3), name::String="", color = RGBA(0.75, 0.75, 0.75)
        )
        T = promote_type(eltype.((r, h, m, xoffset, qoffset))...)
        
        J = m * diagm([1.0; 1.0; 1.0])

        return Body(m, J; name=name, shape=new{T}(xoffset, qoffset, [r; h], scale, color))
    end
end

"""
$(TYPEDEF)

A `Sphere`.
# Important attributes
* `r`: The sphere radius.  
"""
mutable struct Sphere{T} <: Shape{T}
    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    r::T
    scale::SVector{3,T}
    color::RGBA

    function Sphere(r::Real;
            xoffset::AbstractVector = szeros(3), qoffset::UnitQuaternion = one(UnitQuaternion),
            scale::AbstractVector = sones(3), color = RGBA(0.75, 0.75, 0.75)
        )
        T = promote_type(eltype.((r, xoffset, qoffset))...)

        new{T}(xoffset, qoffset, r, scale, color)
    end

    function Sphere(r::Real, m::Real;
            xoffset::AbstractVector = szeros(3), qoffset::UnitQuaternion = one(UnitQuaternion),
            scale::AbstractVector = sones(3), name::String="", color = RGBA(0.75, 0.75, 0.75)
        )
        T = promote_type(eltype.((r, m, xoffset, qoffset))...)
        J = 2 / 5 * m * diagm([r^2 for i = 1:3])

        return Body(m, J; name=name, shape=new{T}(xoffset, qoffset, r, scale, color))
    end
end

"""
$(TYPEDEF)

A square `Pyramid` with a base in the x-y-plane.
* `wh`: The pyramid width and height (vector). 
"""
mutable struct Pyramid{T} <: Shape{T}
    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    wh::SVector{2,T}
    scale::SVector{3,T}
    color::RGBA

    # Pyramid points in the z direction, Center of mass at 1/4 h
    function Pyramid(w::Real, h::Real;
            xoffset::AbstractVector = szeros(3), qoffset::UnitQuaternion = one(UnitQuaternion),
            scale::AbstractVector = sones(3), color = RGBA(0.75, 0.75, 0.75)
        )
        T = promote_type(eltype.((w, h, xoffset, qoffset))...)

        new{T}(xoffset, qoffset, [w;h], scale, color)
    end

    function Pyramid(w::Real, h::Real, m::Real;
            xoffset::AbstractVector = szeros(3), qoffset::UnitQuaternion = one(UnitQuaternion),
            scale::AbstractVector = sones(3), name::String="", color = RGBA(0.75, 0.75, 0.75)
        )
        T = promote_type(eltype.((w, h, m, xoffset, qoffset))...)
        J = 1/80 * m * diagm([4*w^2+3*h^2;4*w^2+3*h^2;8*w^2])

        return Body(m, J; name=name, shape=new{T}(xoffset, qoffset, [w;h], scale, color))
    end
end

# function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, shape::Shape{T}) where {T}
#     summary(io, shape)
#     println(io,"")
#     println(io," xoffset: "*string(shape.xoffset))
#     println(io," qoffset: "*string(shape.qoffset))
#     println(io," color:   "*string(shape.color))
#     println(io," scale:   "*string(shape.scale))
# end

# function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, shape::EmptyShape{T}) where {T}
#     summary(io, shape)
#     println(io,"")
# end
