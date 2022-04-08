"""
    SoftHalfSpaceCollision

    collision between a soft contact and a flat surface

    collider: soft object
    contact_tangent: mapping from world frame to surface tangent frame
    contact_normal: inverse/complement of contact_tangent
    collider_origin: position of contact on Body relative to center of mass
"""
mutable struct SoftHalfSpaceCollision{T,O,I,OI,N} <: SoftCollision{T,O,I,OI,N}
    collider::SoftCollider{T,N}
    contact_normal::Adjoint{T,SVector{I,T}}
    collider_origin::SVector{I,T}
    contact_tangent::SMatrix{O,I,T,OI}
end

function inside(collision::SoftHalfSpaceCollision{T}, p, xc, qc) where T
    origin = szeros(T,3)
    normal = collision.contact_normal
    c = normal * (p - origin)
    return c <= 0.0
end

# normal projection (from child to parent)
function contact_normal(collision::SoftHalfSpaceCollision, p, xc, qc)
    return collision.contact_normal[1,:]
end

parent_origin(collision::SoftHalfSpaceCollision) = collision.collider_origin
child_origin(collision::SoftHalfSpaceCollision{T}) where T = szeros(T,3)
