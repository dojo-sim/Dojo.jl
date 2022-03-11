abstract type Collision{T,O,I,OI} end 

"""
    SphereFloorCollision 

    collision between a spherical contact and a flat surface

    contact_tangent: mapping from world frame to surface tangent frame 
    contact_normal: inverse/complement of contact_tangent
    contact_point: position of contact on Body relative to center of mass 
    contact_radius: radius of contact
"""
mutable struct SphereFloorCollision{T,O,I,OI} <: Collision{T,O,I,OI}
    contact_tangent::SMatrix{O,I,T,OI}
    contact_normal::Adjoint{T,SVector{I,T}}
    contact_point::SVector{I,T}
    contact_radius::T
end 

struct SphereSphereCollision end

struct CapsuleCapsuleCollision end