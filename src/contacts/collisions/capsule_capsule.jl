"""
    CapsuleCapsuleCollision 

    collision between two capsules

    origin_parent: radius of contact for parent body 
    origin_child: radius of contact for child body
    orientation_parent: radius of contact for parent body 
    orientation_child: radius of contact for child body
    radius_parent: radius of contact for parent body 
    radius_child: radius of contact for child body
    height_parent: radius of contact for parent body 
    height_child: radius of contact for child body
"""
mutable struct CapsuleCapsuleCollision{T,O,I,OI} <: Collision{T,O,I,OI}
    origin_parent::SVector{I,T}
    origin_child::SVector{I,T}
    orientation_parent::Quaternion{T}
    orientation_child::Quaternion{T}
    radius_parent::T
    radius_child::T
    height_parent::T
    height_child::T
end 

# distance
function distance(collision::CapsuleCapsuleCollision, xp, qp, xc, qc)
    # contact points
    pp = contact_point(:parent, collision, xp, qp, xc, qc)
    pc = contact_point(:child,  collision, xp, qp, xc, qc)

    # difference
    d = pp - pc
     
    return norm(pp - pc, 2)
end

function ∂distance∂x(gradient::Symbol, collision::CapsuleCapsuleCollision, xp, qp, xc, qc)
    # contact points
    pp = contact_point(:parent, collision, xp, qp, xc, qc)
    pc = contact_point(:child,  collision, xp, qp, xc, qc)

    # difference
    d = pp - pc

    X = ∂norm∂x(d) * (∂contact_point∂x(:parent, jacobian, collision, xp, qp, xc, qc) - ∂contact_point∂x(:child, jacobian, collision, xp, qp, xc, qc))

    return X
end

function ∂distance∂q(gradient::Symbol, collision::CapsuleCapsuleCollision, xp, qp, xc, qc)
    # contact points
    pp = contact_point(:parent, collision, xp, qp, xc, qc)
    pc = contact_point(:child,  collision, xp, qp, xc, qc)

    # difference
    d = pp - pc

    X = ∂norm∂x(d) * (∂contact_point∂q(:parent, jacobian, collision, xp, qp, xc, qc) - ∂contact_point∂q(:child, jacobian, collision, xp, qp, xc, qc))

    return X
end

# contact point in world frame
function contact_point(relative::Symbol, collision::CapsuleCapsuleCollision, xp, qp, xc, qc) 
    # call mehrotra 

    if relative == :parent 
        return collision.ip.z[0 .+ (1:3)]
    elseif relative == :child
        return collision.ip.z[3 .+ (1:3)]
    end
end

function ∂contact_point∂x(relative::Symbol, jacobian::Symbol, collision::CapsuleCapsuleCollision, xp, qp, xc, qc)
    # call mehrotra 

    if relative == :parent 
        if jacobian == :parent
            X = collision.ip.δz[]
        elseif jacobian == :child 
            X = collision.ip.δz[]
        end
    elseif relative == :child
        if jacobian == :parent
            X = collision.ip.δz[]
        elseif jacobian == :child 
            X = collision.ip.δz[]
        end
    end

    return X
end

function ∂contact_point∂q(relative::Symbol, jacobian::Symbol, collision::CapsuleCapsuleCollision, xp, qp, xc, qc)
# call mehrotra 

    if relative == :parent 
        if jacobian == :parent
            X = collision.ip.δz[]
        elseif jacobian == :child 
            X = collision.ip.δz[]
        end
    elseif relative == :child
        if jacobian == :parent
            X = collision.ip.δz[]
        elseif jacobian == :child 
            X = collision.ip.δz[]
        end
    end
    
    return X
end

