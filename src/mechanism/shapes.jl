abstract type Shape{T} end

struct EmptyShape{T} <: Shape{T}
    EmptyShape() = new{Float64}()
end

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

        mass_cylinder = π * h * r^2.0
        mass_hemisphere = π * 2.0 / 3.0 * r^3.0 
        mass_total = mass_cylinder + 2 * mass_hemisphere
        Ixx_cylinder = mass_cylinder * (1.0 / 12.0 * h^2.0 + 1.0 / 4.0 * r^2.0)
        Ixx_hemisphere = 83.0 / 320 * mass_hemisphere * r^2
        d = (3.0 / 8.0 * r + 0.5 * h)
        Ixx = Ixx_cylinder + 2.0 * (Ixx_hemisphere + mass_hemisphere * d^2.0)
        Izz = 0.5 * mass_cylinder * r^2.0 + mass_hemisphere * (2.0 * r^2.0) / 5.0 

        J = m * diagm([Ixx; Ixx; Izz])

        return Body(m, J; name=name, shape=new{T}(xoffset, qoffset, [r; h], scale, color))
    end
end

mutable struct Shapes14{T} <: Shape{T}
    shape::Vector 
    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}
    scale::SVector{3,T}
    color::RGBA

    function Shapes14(shapes::Vector; 
        xoffset::AbstractVector = szeros(3), qoffset::UnitQuaternion = one(UnitQuaternion),
        scale::AbstractVector = sones(3), name::String="", color = RGBA(0.75, 0.75, 0.75)) where {T}
        new{T}(shapes, xoffset, qoffset, scale, color)
    end

    function Shapes14(shapes::Vector, m::T, J; xoffset::AbstractVector = szeros(3), qoffset::UnitQuaternion = one(UnitQuaternion),
        scale::AbstractVector = sones(3), name::String="", color = RGBA(0.75, 0.75, 0.75)) where T
        Body(m, J; name=name, shape=new{T}(shapes, xoffset, qoffset, scale, color))
    end
end

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
