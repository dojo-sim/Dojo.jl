"""
    SphereFlatCollision 

    collision between a spherical contact and a flat surface

    contact_tangent: mapping from world frame to surface tangent frame 
    contact_normal: inverse/complement of contact_tangent
    contact_origin: position of contact on Body relative to center of mass 
    contact_radius: radius of contact
"""
mutable struct SphereFlatCollision{T,O,I,OI} <: Collision{T,O,I,OI}
    contact_tangent::SMatrix{O,I,T,OI}
    contact_normal::Adjoint{T,SVector{I,T}}
    contact_origin::SVector{I,T}
    contact_radius::T
end 

# distance
function distance(collision::SphereFlatCollision, xp, qp, xc, qc)
    collision.contact_normal * (xp + vector_rotate(collision.contact_origin, qp)) - collision.contact_radius
end

function ∂distance∂x(relative::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc)
    if relative == :parent
        return collision.contact_normal
    elseif relative == :child 
        return szeros(eltype(xc), 1, 3) 
    end
end

function ∂distance∂q(relative::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc)
    if relative == :parent
        return collision.contact_normal * ∂vector_rotate∂q(collision.contact_origin, qp)
    elseif relative == :child 
        return szeros(eltype(qc), 1, 4)
    end
end

# contact point in world frame
function contact_point(relative::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc) 
    if relative == :parent
        return xp + vector_rotate(model.collision.contact_origin, qp) - model.collision.contact_normal' * model.collision.contact_radius
    elseif relative == :child 
        projector = model.collision.contact_tangent' * model.collision.contact_tangent 
        return projector * (xp + vector_rotate(model.collision.contact_origin, qp))
    end
end

function ∂contact_point∂x(relative::Symbol, jacobian::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc)
    if relative == :parent 
        if jacobian == :parent 
            return 1.0 * I(3)
        elseif jacobian == :child 
            return szeros(eltype(xp), 3, 3)
        end
    elseif relative == :child 
        if jacobian == :parent 
            projector = model.collision.contact_tangent' * model.collision.contact_tangent 
            return projector
        elseif jacobian == :child 
            return szeros(eltype(xp), 3, 3)
        end
    end
end

function ∂contact_point∂q(relative::Symbol, jacobian::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc)
    if relative == :parent 
        if jacobian == :parent 
            return ∂vector_rotate∂q(model.collision.contact_origin, qp)
        end
    elseif relative == :child 
        if jacobian == :parent 
            projector = model.collision.contact_tangent' * model.collision.contact_tangent 
            return projector * ∂vector_rotate∂q(model.collision.contact_origin, qp)
        elseif jacobian == :child 
            return szeros(eltype(qp), 3, 4)
        end
    end
end
    
# normal projection (from child to parent)
function contact_normal(relative::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc)
    if relative == :parent
        return collision.contact_normal
    elseif relative == :child 
        return -1.0 * collision.contact_normal
    end
end

function ∂contact_normal∂x(relative::Symbol, jacobian::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc, λ)
    contact_normal = collision.contact_normal
    no, ni = size(contact_normal)
    return szeros(eltype(contact_normal), no, 3)
end

function ∂contact_normal∂q(relative::Symbol, jacobian::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc, λ)
    contact_normal = collision.contact_normal
    no, ni = size(contact_normal)
    return szeros(eltype(contact_normal), no, 4)
end

# tangent projection
function contact_tangent(relative::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc)
    if relative == :parent 
        return collision.contact_tangent
    elseif relative == :child 
        return -collision.contact_tangent 
    end
end

function ∂contact_tangent∂x(relative::Symbol, jacobian::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc, λ)
    contact_tangent = collision.contact_tangent
    no, ni = size(contact_tangent)
    return szeros(eltype(contact_tangent), no, 3)
end

function ∂contact_tangent∂q(relative::Symbol, jacobian::Symbol, collision::SphereFlatCollision, xp, qp, xc, qc, λ)
    contact_tangent = collision.contact_tangent
    no, ni = size(contact_tangent)
    return szeros(eltype(contact_tangent), no, 4)
end

