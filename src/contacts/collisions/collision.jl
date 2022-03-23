"""
    Collision

    abstract type defining interaction between two bodies
"""
abstract type Collision{T,O,I,OI} end

# contact point origin
function contact_point_origin(x, q, k)
    x + vector_rotate(k, q)
end

∂contact_point_origin∂x(x, q, k) = 1.0 * I(3)
∂contact_point_origin∂q(x, q, k) = ∂vector_rotate∂q(k, q)
∂contact_point_origin∂k(x, q, k) = rotation_matrix(q)

"""
    contact_normal(collision, xp, qp, xc, qc)

    the contact normal (from child to parent) between two contact points

    collision: Collision
    xp: parent body position
    qp: parent body orientation
    xc: child body position
    qc: child body orientation
"""
function contact_normal(collision::Collision, xp, qp, xc, qc)
    # contact points
    cop = contact_point(:parent, collision, xp, qp, xc, qc)
    coc = contact_point(:child,  collision, xp, qp, xc, qc)

    # unnormalized direction
    dir = cop - coc

    # distance
    dis = distance(collision, xp, qp, xc, qc)

    # normalized direction
    if dis >= 0.0
        return normalize(dir)'
    else
        return -1.0 * normalize(dir)'
    end
end

function ∂contact_normal_transpose∂x(jacobian::Symbol, collision::Collision, xp, qp, xc, qc)
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
    if dis >= 0.0
        return X
    else
        return -1.0 * X
    end
end

function ∂contact_normal_transpose∂q(jacobian::Symbol, collision::Collision, xp, qp, xc, qc)
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
    if dis >= 0.0
        return Q
    else
        return -1.0 * Q
    end
end

function contact_tangent_one(collision::Collision, xp, qp, xc, qc)
    # normal
    n = contact_normal(collision, xp, qp, xc, qc)

    # tangents
    w = [1.0; 0.0; 0.0]
    t1 = skew(w) * n' # tangent
    if !(norm(t1) > 1.0e-6)
        # @warn "edge case!"
        w = [0.0; 1.0; 0.0]
        t1 = skew(w) * n' # tangent
    end

    return t1'
end

function contact_tangent_two(collision::Collision, xp, qp, xc, qc)
    # normal
    n = contact_normal(collision, xp, qp, xc, qc)

    # tangent one
    t1 = contact_tangent_one(collision, xp, qp, xc, qc)

    # tangent two
    t2 = skew(t1') * n'

    return t2'
end

"""
    contact_tangent(collision, xp, qp, xc, qc)

    contact tangents between two contact points

    collision: Collision
    xp: parent body position
    qp: parent body orientation
    xc: child body position
    qc: child body orientation
"""
function contact_tangent(collision::Collision, xp, qp, xc, qc)
   [
       contact_tangent_one(collision, xp, qp, xc, qc);
       contact_tangent_two(collision, xp, qp, xc, qc);
   ]
end

function ∂contact_tangent_one_transpose∂x(jacobian::Symbol, collision::Collision, xp, qp, xc, qc)
    # normal
    n = contact_normal(collision, xp, qp, xc, qc)

    # tangents
    w = [1.0; 0.0; 0.0] # candidate
    t1 = skew(w) * n' # tangent
    if !(norm(t1) > 1.0e-5)
        w = [0.0; 1.0; 0.0] # candidate
        t1 = skew(w) * n' # tangent
    end

    ∂t1∂x = skew(w) * ∂contact_normal_transpose∂x(jacobian, collision, xp, qp, xc, qc)

    # if jacobian == :parent
    #     FD = FiniteDiff.finite_difference_jacobian(x -> contact_tangent(collision, x, qp, xc, qc)[1, :]', xp)
    # elseif jacobian == :child
    #     FD = FiniteDiff.finite_difference_jacobian(x -> contact_tangent(collision, xp, qp, x, qc)[1, :]', xc)
    # end

    # @assert norm(FD - ∂t1∂x, Inf) < 1.0e-2

    return ∂t1∂x
end

function ∂contact_tangent_two_transpose∂x(jacobian::Symbol, collision::Collision, xp, qp, xc, qc)
    # normal
    n = contact_normal(collision, xp, qp, xc, qc)
    t1 = contact_tangent_one(collision, xp, qp, xc, qc)

    ∂t1∂x = ∂contact_tangent_one_transpose∂x(jacobian, collision, xp, qp, xc, qc)
    ∂t2∂x = skew(t1') * ∂contact_normal_transpose∂x(jacobian, collision, xp, qp, xc, qc) + ∂skew∂p(n') * ∂t1∂x

    # if jacobian == :parent
    #     FD = FiniteDiff.finite_difference_jacobian(x -> contact_tangent(collision, x, qp, xc, qc)[2, :]', xp)
    # elseif jacobian == :child
    #     FD = FiniteDiff.finite_difference_jacobian(x -> contact_tangent(collision, xp, qp, x, qc)[2, :]', xc)
    # end

    # @assert norm(FD - ∂t2∂x, Inf) < 1.0e-2

    return ∂t2∂x
end

function ∂contact_tangent_one_transpose∂q(jacobian::Symbol, collision::Collision, xp, qp, xc, qc)
    # normal
    n = contact_normal(collision, xp, qp, xc, qc)

    # tangents
    w = [1.0; 0.0; 0.0] # candidate
    t1 = skew(w) * n' # tangent
    if !(norm(t1) > 1.0e-5)
        w = [0.0; 1.0; 0.0] # candidate
        t1 = skew(w) * n' # tangent
    end

    ∂t1∂q = skew(t1) * ∂contact_normal_transpose∂q(jacobian, collision, xp, qp, xc, qc)

    # if jacobian == :parent
    #     FD = FiniteDiff.finite_difference_jacobian(q -> contact_tangent(collision, xp, Quaternion(q..., false), xc, qc)[1, :]', vector(qp))
    # elseif jacobian == :child
    #     FD = FiniteDiff.finite_difference_jacobian(q -> contact_tangent(collision, xp, qp, xc, Quaternion(q..., false))[1, :]', vector(qc))
    # end

    # return FD
    # @assert norm(FD - ∂t1∂q, Inf) < 1.0e-2

    return ∂t1∂q
end

function ∂contact_tangent_two_transpose∂q(jacobian::Symbol, collision::Collision, xp, qp, xc, qc)
    # normal
    n = contact_normal(collision, xp, qp, xc, qc)

    # tangents
    w = [1.0; 0.0; 0.0] # candidate
    t1 = skew(w) * n' # tangent
    if !(norm(t1) > 1.0e-5)
        w = [0.0; 1.0; 0.0] # candidate
        t1 = skew(w) * n' # tangent
    end
    t2 = skew(t1) * n' # tangent

    ∂t1∂q = skew(w) * ∂contact_normal_transpose∂q(jacobian, collision, xp, qp, xc, qc)

    ∂t2∂q = skew(t1) * ∂contact_normal_transpose∂q(jacobian, collision, xp, qp, xc, qc) + ∂skew∂p(n') * ∂t1∂q

    # if jacobian == :parent
    #     FD = FiniteDiff.finite_difference_jacobian(q -> contact_tangent(collision, xp, Quaternion(q..., false), xc, qc)[2, :]', vector(qp))
    # elseif jacobian == :child
    #     FD = FiniteDiff.finite_difference_jacobian(q -> contact_tangent(collision, xp, qp, xc, Quaternion(q..., false))[2, :]', vector(qc))
    # end

    # @assert norm(FD - ∂t2∂q, Inf) < 1.0e-2

    return ∂t2∂q
end
