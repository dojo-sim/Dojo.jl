"""
    SphereFlatCollision 

    collision between a spherical contact and a flat surface

    contact_tangent: mapping from world frame to surface tangent frame 
    contact_normal: inverse/complement of contact_tangent
    contact_point: position of contact on Body relative to center of mass 
    contact_radius: radius of contact
"""
mutable struct SphereFlatCollision{T,O,I,OI} <: Collision{T,O,I,OI}
    contact_tangent::SMatrix{O,I,T,OI}
    contact_normal::SVector{I,T}
    contact_point::SVector{I,T}
    contact_radius::T
end 

# distance
function distance(collision::SphereFlatCollision, xp, qp, xc, qc)
    collision.contact_normal' * (xp + vector_rotate(collision.contact_point, qp)) - collision.contact_radius
end

function ∂distance∂x(relative::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc)
    collision.contact_normal'
end

function ∂distance∂q(relative::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc)
    collision.contact_normal' * ∂vector_rotate∂q(collision.contact_point, qp)
end

# contact point in body frame
function contact_point(collision::SphereFlatCollision, xp, qp, xc, qc) 
    return collision.contact_point 
end

function ∂contact_point∂x(relative::Symbol, jacobian::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc)
    szeros(eltype(collision.contact_point), 3, 3)
end

function ∂contact_point∂q(relative::Symbol, jacobian::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc)
    szeros(eltype(collision.contact_point), 3, 4)
end

# contact radius 
function contact_radius(relative::Symbol, collision::SphereFlatCollision)
    if relative == :parent 
        return collision.contact_radius
    elseif relative == :child 
        return 0.0
    end
end
    
# normal projection 
function contact_normal(collision::SphereFlatCollision, xp, qp, xc, qc)
    collision.contact_normal
end

function ∂contact_normal∂x(jacobian::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc, Δv)
    contact_normal = collision.contact_normal
    no, ni = size(contact_normal)
    return szeros(eltype(contact_normal), no, 3)
end

function ∂contact_normal∂q(jacobian::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc, Δv)
    contact_normal = collision.contact_normal
    no, ni = size(contact_normal)
    return szeros(eltype(contact_normal), no, 4)
end

# tangent projection
function contact_tangent(collision::SphereFlatCollision, xp, qp, xc, qc)
    collision.contact_tangent
end

function ∂contact_tangent∂x(jacobian::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc, Δv)
    contact_tangent = collision.contact_tangent
    no, ni = size(contact_tangent)
    return szeros(eltype(contact_tangent), no, 3)
end

function ∂contact_tangent∂q(jacobian::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc, Δv)
    contact_tangent = collision.contact_tangent
    no, ni = size(contact_tangent)
    return szeros(eltype(contact_tangent), no, 4)
end

