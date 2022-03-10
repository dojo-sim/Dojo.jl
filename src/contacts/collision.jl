abstract type Collision{T,O,I,OI} end 

"""
    SphereFloorCollision 

    collision between a spherical contact and a flat surface

    surface_projector: mapping from world frame to surface tangent frame 
    surface_normal_projector: inverse/complement of surface_projector
    contact_point: position of contact on Body relative to center of mass 
    contact_radius: radius of contact
"""
mutable struct SphereFloorCollision{T,O,I,OI} <: Collision{T,O,I,OI}
    surface_projector::SMatrix{O,I,T,OI}
    surface_normal_projector::Adjoint{T,SVector{I,T}}
    contact_point::SVector{I,T}
    contact_radius::T
end 

struct SphereSphereCollision end

struct CapsuleCapsuleCollision end