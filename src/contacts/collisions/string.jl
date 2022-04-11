"""
    StringCollision 

    string collision between contact points

    origin_parent: position of contact on parent body relative to center of mass 
    origin_child: position of contact on parent body relative to center of mass 
    length: maximum distance between contact point
"""
mutable struct StringCollision{T,O,I,OI} <: Collision{T,O,I,OI}
    origin_parent::SVector{I,T}
    origin_child::SVector{I,T}
    length::T
end 

# distance
function distance(collision::StringCollision, xp, qp, xc, qc)
    # contact origin points
    cop = contact_point_origin(xp, qp, collision.origin_parent) 
    coc = contact_point_origin(xc, qc, collision.origin_child)

    # distance between contact origins
    d = norm(coc - cop, 2)

    # minimum distance between spheres
    return collision.length - d
end

function ∂distance∂x(gradient::Symbol, collision::StringCollision, xp, qp, xc, qc)
    # contact origin points
    cop = contact_point_origin(xp, qp, collision.origin_parent) 
    coc = contact_point_origin(xc, qc, collision.origin_child)

    # distance between contact origins
    # d = norm(cop - coc, 2)
    ∂norm∂d = ∂norm∂x(coc - cop)

    if gradient == :parent
        D = 1.0 * ∂norm∂d *  -1.0 * ∂contact_point_origin∂x(xp, qp, collision.origin_parent)
    elseif gradient == :child 
        D = 1.0 * ∂norm∂d * 1.0 * ∂contact_point_origin∂x(xc, qc, collision.origin_child)
    end

    if gradient == :parent 
        FD = FiniteDiff.finite_difference_jacobian(x -> distance(collision, x, qp, xc, qc), xp)
    elseif gradient == :child 
        FD = FiniteDiff.finite_difference_jacobian(x -> distance(collision, xp, qp, x, qc), xc)
    end

    @assert norm(-D - FD, Inf) < 1.0e-5

    return -D
end

function ∂distance∂q(gradient::Symbol, collision::StringCollision, xp, qp, xc, qc)
    # contact origin points
    cop = contact_point_origin(xp, qp, collision.origin_parent) 
    coc = contact_point_origin(xc, qc, collision.origin_child)

    # distance between contact origins
    # d = norm(cop - coc, 2)
    ∂norm∂d = ∂norm∂x(coc - cop)
    
    if gradient == :parent
        D = 1.0 * ∂norm∂d *  -1.0 * ∂contact_point_origin∂q(xp, qp, collision.origin_parent)
    elseif gradient == :child 
        D = 1.0 * ∂norm∂d * 1.0 * ∂contact_point_origin∂q(xc, qc, collision.origin_child)
    end

    if gradient == :parent 
        FD = FiniteDiff.finite_difference_jacobian(q -> distance(collision, xp, Quaternion(q..., false), xc, qc), vector(qp))
    elseif gradient == :child 
        FD = FiniteDiff.finite_difference_jacobian(q -> distance(collision, xp, qp, xc, Quaternion(q..., false)), vector(qc))
    end

    @assert norm(-D - FD, Inf) < 1.0e-5

    return -D
end

# contact point in world frame
function contact_point(relative::Symbol, collision::StringCollision, xp, qp, xc, qc) 
    # contact origin points
    cop = contact_point_origin(xp, qp, collision.origin_parent) 
    coc = contact_point_origin(xc, qc, collision.origin_child)

    # direction of minimum distance (child to parent)
    # d = cop - coc 
    # dir = normalize(d)

    # contact point
    if relative == :parent
        return cop
    elseif relative == :child 
        return coc
    end
end

function ∂contact_point∂x(relative::Symbol, jacobian::Symbol, collision::StringCollision, xp, qp, xc, qc)
    # contact origin points
    cop = contact_point_origin(xp, qp, collision.origin_parent) 
    coc = contact_point_origin(xc, qc, collision.origin_child)

    # direction of minimum distance (child to parent)
    # d = cop - coc 
    # dir = normalize(d)

    if relative == :parent 
        # cop
        if jacobian == :parent 
            ∂c∂x = ∂contact_point_origin∂x(xp, qp, collision.origin_parent)
            X = ∂c∂x 
        elseif jacobian == :child 
            X = szeros(eltype(xc), 3, 3)
        end
    elseif relative == :child 
        # coc
        if jacobian == :parent 
            X = szeros(eltype(xp), 3, 3)
        elseif jacobian == :child 
            ∂c∂x = ∂contact_point_origin∂x(xc, qc, collision.origin_child)
            X = ∂c∂x 
        end
    end

    if jacobian == :parent
        FD =  FiniteDiff.finite_difference_jacobian(x -> contact_point(relative, collision, x, qp, xc, qc), xp)
    elseif jacobian == :child 
        FD = FiniteDiff.finite_difference_jacobian(x -> contact_point(relative, collision, xp, qp, x, qc), xc)
    end

    @assert norm(X - FD, Inf) < 1.0e-5

    return X
end

function ∂contact_point∂q(relative::Symbol, jacobian::Symbol, collision::StringCollision, xp, qp, xc, qc)
    # contact origin points
    cop = contact_point_origin(xp, qp, collision.origin_parent) 
    coc = contact_point_origin(xc, qc, collision.origin_child)

    # direction of minimum distance (child to parent)
    # d = cop - coc 
    # dir = normalize(d)

    if relative == :parent 
        # cop
        if jacobian == :parent 
            ∂c∂q = ∂contact_point_origin∂q(xp, qp, collision.origin_parent)
            Q = ∂c∂q 
        elseif jacobian == :child 
            Q = szeros(eltype(xc), 3, 4)
        end
    elseif relative == :child 
        # coc
        if jacobian == :parent 
            Q = szeros(eltype(xp), 3, 4)
        elseif jacobian == :child 
            ∂c∂q = ∂contact_point_origin∂q(xc, qc, collision.origin_child)
            Q = ∂c∂q 
        end
    end

    if jacobian == :parent
        FD = FiniteDiff.finite_difference_jacobian(q -> contact_point(relative, collision, xp, Quaternion(q..., false), xc, qc), vector(qp))
    elseif jacobian == :child 
        FD = FiniteDiff.finite_difference_jacobian(q -> contact_point(relative, collision, xp, qp, xc, Quaternion(q..., false)), vector(qc))
    end

    @assert norm(Q - FD, Inf) < 1.0e-5

    return Q
end

function contact_normal(collision::StringCollision, xp, qp, xc, qc)
    # contact points
    cop = contact_point(:parent, collision, xp, qp, xc, qc)
    coc = contact_point(:child,  collision, xp, qp, xc, qc)

    # unnormalized direction
    dir = cop - coc

    # distance
    dis = distance(collision, xp, qp, xc, qc)

    # normalized direction
    # if dis >= 0.0
    return -normalize(dir)'
    # else
    #     return normalize(dir)'
    # end
end

function ∂contact_normal_transpose∂x(jacobian::Symbol, collision::StringCollision, xp, qp, xc, qc)
    # contact origin points
    cop = contact_point(:parent, collision, xp, qp, xc, qc)
    coc = contact_point(:child,  collision, xp, qp, xc, qc)

    # unnormalized direction
    dir = cop - coc

    # Jacobians
    X = ∂normalize∂x(dir) * (∂contact_point∂x(:parent, jacobian, collision, xp, qp, xc, qc) - ∂contact_point∂x(:child, jacobian, collision, xp, qp, xc, qc))

    # distance
    dis = distance(collision, xp, qp, xc, qc)

    # normalized direction
    # if dis >= 0.0
    return -X
    # else
    #     return 1.0 * X
    # end
end

function ∂contact_normal_transpose∂q(jacobian::Symbol, collision::StringCollision, xp, qp, xc, qc)
    # contact origin points
    cop = contact_point(:parent, collision, xp, qp, xc, qc)
    coc = contact_point(:child,  collision, xp, qp, xc, qc)

    # unnormalized direction
    dir = cop - coc

    Q = ∂normalize∂x(dir) * (∂contact_point∂q(:parent, jacobian, collision, xp, qp, xc, qc) - ∂contact_point∂q(:child, jacobian, collision, xp, qp, xc, qc))

    # Jacobians
    # if jacobian == :parent
    #     FD = FiniteDiff.finite_difference_jacobian(q -> contact_normal(collision, xp, Quaternion(q..., false), xc, qc)', vector(qp))
    # elseif jacobian == :child
    #     FD = FiniteDiff.finite_difference_jacobian(q -> contact_normal(collision, xp, qp, xc, Quaternion(q..., false))', vector(qc))
    # end

    # distance
    dis = distance(collision, xp, qp, xc, qc)

    # @assert norm((dis >= 0.0 ? 1.0 : -1.0) * Q - FD, Inf) < 1.0e-2

    # normalized direction
    # if dis >= 0.0
    return -Q
    # else
    #     return Q
    # end
end