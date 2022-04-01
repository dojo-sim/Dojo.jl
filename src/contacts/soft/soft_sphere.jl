"""
    SoftSphereCollision

    collision between a soft contact and a sphere

    collider: soft object
    contact_tangent: mapping from world frame to surface tangent frame
    collider_origin: position of contact on collider relative to its center of mass
    sphere_origin: position of contact on sphere relative to its center of mass
    contact_radius: radius of contact
"""
mutable struct SoftSphereCollision{T,O,I,OI,N} <: SoftCollision{T,O,I,OI,N}
    collider::SoftCollider{T,N}
    contact_tangent::SMatrix{O,I,T,OI}
    collider_origin::SVector{I,T}
    sphere_origin::SVector{I,T}
    contact_radius::T
end

function inside(collision::SoftSphereCollision{T}, p, xc, qc) where T
    # p = barycenter of the overlap expressed in the world frame
    # pc = center of the sphere expressed in the world frame
    pc = xc + Dojo.vector_rotate(collision.sphere_origin, qc)
    c = norm(p - pc)
    return c <= collision.contact_radius
end

# normal projection (from child to parent)
function contact_normal(collision::SoftSphereCollision, p, xc, qc)
    # p = barycenter of the overlap expressed in the world frame
    # pc = center of the sphere expressed in the world frame
    pc = xc + Dojo.vector_rotate(collision.sphere_origin, qc)
    normal = p - pc
    return normal ./ (1e-20 + norm(normal))
end
