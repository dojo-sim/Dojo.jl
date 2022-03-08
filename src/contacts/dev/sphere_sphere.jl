using LinearAlgebra 
using MeshCat
using GeometryBasics
using Colors
using CoordinateTransformations

vis = Visualizer()
open(vis)

struct SphereContact{T}
    position::Vector{T}
    radius::T
end

p1 = [-0.11, 0.0, 0.0]
p2 = [0.11, 0.0, 0.0]
s1 = SphereContact(p1, 0.1)
s2 = SphereContact(p2, 0.1)
m1 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), s1.radius)
m2 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), s2.radius)
setobject!(vis[:s1], m1, MeshPhongMaterial(color=Colors.RGBA(1.0, 0.0, 0.0, 1.0)))
setobject!(vis[:s2], m2, MeshPhongMaterial(color=Colors.RGBA(1.0, 0.0, 0.0, 1.0)))
settransform!(vis[:s1], Translation(s1.position))
settransform!(vis[:s2], Translation(s2.position))

function sdf_sphere_sphere(s1::SphereContact, s2::SphereContact)
    position_error = s2.position - s1.position 
    total_radius = s1.radius + s2.radius
    norm(position_error) - total_radius # should be >= 0
end

sdf_sphere_sphere(s1, s2)

function sdf_sq_sphere_sphere(s1::SphereContact, s2::SphereContact)
    position_error = s2.position - s1.position 
    total_radius = s1.radius + s2.radius
    dot(position_error, position_error) - total_radius^2 # should be >= 0
end

sdf_sq_sphere_sphere(s1, s2)

function direction_sphere_sphere(s1::SphereContact, s2::SphereContact)
    # from sphere 1 to sphere 2
    position_error = s2.position - s1.position 
    magnitude = norm(position_error)

    if magnitude > 0.0
        direction = position_error ./ magnitude
    else 
        direction = [1.0, 0.0, 0.0]
    end

    return direction
end

direction_sphere_sphere(s1, s2)

function closest_points_sphere_sphere(s1::SphereContact, s2::SphereContact)
    # direction 
    direction = direction_sphere_sphere(s1::SphereContact, s2::SphereContact)

    # closest point on sphere 1
    p1 = s1.position + direction * s1.radius 

    # closest point on sphere 2
    p2 = s2.position - direction * s2.radius 

    return p1, p2
end

c1, c2 = closest_points_sphere_sphere(s1, s2)

cs1 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.01)
cs2 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.01)
setobject!(vis[:cs1], cs1, MeshPhongMaterial(color=Colors.RGBA(0.0, 0.0, 0.0, 1.0)))
setobject!(vis[:cs2], cs2, MeshPhongMaterial(color=Colors.RGBA(0.0, 0.0, 0.0, 1.0)))
settransform!(vis[:cs1], Translation(c1))
settransform!(vis[:cs2], Translation(c2))
