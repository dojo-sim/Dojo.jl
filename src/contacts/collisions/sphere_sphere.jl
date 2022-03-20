"""
    SphereSphereCollision 

    collision between two spheres

    # contact_tangent: directions tangent to contact normal 
    # contact_normal: direction of minimum distance between to contacts (child to parent)
    contact_origin_parent: position of contact on parent body relative to center of mass 
    contact_origin_child: position of contact on parent body relative to center of mass 
    contact_radius_parent: radius of contact for parent body 
    contact_radius_child: radius of contact for child body
"""
mutable struct SphereSphereCollision{T,O,I,OI} <: Collision{T,O,I,OI}
    contact_origin_parent::SVector{I,T}
    contact_origin_child::SVector{I,T}
    contact_radius_parent::T
    contact_radius_child::T
end 

# contact point origin 
function contact_point_origin(x, q, k) 
    x + vector_rotate(k, q)
end

∂contact_point_origin∂x(x, q, k) = 1.0 * I(3)
∂contact_point_origin∂q(x, q, k) = ∂vector_rotate∂q(k, q)
∂contact_point_origin∂k(x, q, k) = ∂vector_rotate∂p(k, q)

# distance
function distance(collision::SphereSphereCollision, xp, qp, xc, qc)
    # contact origin points
    cop = contact_point_origin(xp, qp, collision.contact_origin_parent) 
    coc = contact_point_origin(xc, qc, collision.contact_origin_child)

    # distance between contact origins
    d = norm(cop - coc, 2)

    # minimum distance between spheres
    return d - (collision.contact_radius_parent + collision.contact_radius_child)
end

function ∂distance∂x(gradient::Symbol, collision::SphereSphereCollision, xp, qp, xc, qc)
    # # contact origin points
    # cop = contact_point_origin(xp, qp, collision.contact_origin_parent) 
    # coc = contact_point_origin(xc, qc, collision.contact_origin_child)

    # # distance between contact origins
    # d = norm(cop - coc, 2)
    # ∂norm∂d = ∂norm∂x(d)

    # if gradient == :parent
    #     return ∂norm∂d *  1.0 * ∂contact_point_origin∂x(xp, qp, collision.contact_origin_parent)
    # elseif gradient == :child 
    #     return ∂norm∂d * -1.0 * ∂contact_point_origin∂x(xc, qc, collision.contact_origin_child)
    # end
    if gradient == :parent 
        return FiniteDiff.finite_difference_jacobian(x -> distance(collision, x, qp, xc, qc), xp)
    elseif gradient == :child 
        return FiniteDiff.finite_difference_jacobian(x -> distance(collision, xp, qp, x, qc), xc)
    end
end

function ∂distance∂q(gradient::Symbol, collision::SphereSphereCollision, xp, qp, xc, qc)
    # contact origin points
    # cop = contact_point_origin(xp, qp, collision.contact_origin_parent) 
    # coc = contact_point_origin(xc, qc, collision.contact_origin_child)

    # # distance between contact origins
    # d = norm(cop - coc, 2)
    # ∂norm∂d = ∂norm∂x(d)
    
    # if gradient == :parent
    #     return ∂norm∂d *  1.0 * ∂contact_point_origin∂q(xp, qp, collision.contact_origin_parent)
    # elseif gradient == :child 
    #     return ∂norm∂d * -1.0 * ∂contact_point_origin∂x(xc, qc, collision.contact_origin_child)
    # end
    if gradient == :parent 
        return FiniteDiff.finite_difference_jacobian(q -> distance(collision, xp, UnitQuaternion(q..., false), xc, qc), vector(qp))
    elseif gradient == :child 
        return FiniteDiff.finite_difference_jacobian(q -> distance(collision, xp, qp, xc, UnitQuaternion(q..., false)), vector(qc))
    end
end

# contact point in world frame
function contact_point(relative::Symbol, collision::SphereSphereCollision, xp, qp, xc, qc) 
    # contact origin points
    cop = contact_point_origin(xp, qp, collision.contact_origin_parent) 
    coc = contact_point_origin(xc, qc, collision.contact_origin_child)

    # direction of minimum distance (child to parent)
    d = cop - coc 
    dir = normalize(d)

    # contact point
    if relative == :parent
        return cop - collision.contact_radius_parent * dir
    elseif relative == :child 
        return coc + collision.contact_radius_child * dir
    end
end

function ∂contact_point∂x(relative::Symbol, jacobian::Symbol, collision::SphereSphereCollision, xp, qp, xc, qc)

    if jacobian == :parent
        return FiniteDiff.finite_difference_jacobian(x -> contact_point(relative, collision, x, qp, xc, qc), xp)
    elseif jacobian == :child 
        return FiniteDiff.finite_difference_jacobian(x -> contact_point(relative, collision, xp, qp, x, qc), xc)
    end
    # # contact origin points
    # cop = contact_point_origin(xp, qp, collision.contact_origin_parent) 
    # coc = contact_point_origin(xc, qc, collision.contact_origin_child)

    # # direction of minimum distance (child to parent)
    # d = cop - coc 
    # dir = normalize(d)

    # if relative == :parent 
    #     # cop - collision.contact_radius_parent * dir
    #     if jacobian == :parent 
    #         ∂c∂x = ∂contact_point_origin∂x(xp, qp, collision.contact_origin_parent)
    #         X = ∂c∂x 
    #         X -= collision.contact_radius_parent * ∂normalize∂x(d) * ∂c∂x
    #         return X
    #     elseif jacobian == :child 
    #         X = -1.0 * collision.contact_radius_parent * ∂normalize∂x(d) * -1.0 * ∂contact_point_origin∂x(xc, qc, collision.contact_origin_child)
    #         return X
    #     end
    # elseif relative == :child 
    #     # coc + collision.contact_radius_child * dir
    #     if jacobian == :parent 
    #         X = 1.0 * collision.contact_radius_child * ∂normalize∂x(d) * ∂contact_point_origin∂x(xp, qp, collision.contact_origin_parent)
    #         return 
    #     elseif jacobian == :child 
    #         ∂c∂x = ∂contact_point_origin∂x(xc, qc, collision.contact_origin_child)
    #         X = ∂c∂x 
    #         X += collision.contact_radius_child * ∂normalize∂x(d) * -1.0 * ∂c∂x
    #         return X
    #     end
    # end
end

function ∂contact_point∂q(relative::Symbol, jacobian::Symbol, collision::SphereSphereCollision, xp, qp, xc, qc)
    if jacobian == :parent
        return FiniteDiff.finite_difference_jacobian(q -> contact_point(relative, collision, xp, UnitQuaternion(q..., false), xc, qc), vector(qp))
    elseif jacobian == :child 
        return FiniteDiff.finite_difference_jacobian(q -> contact_point(relative, collision, xp, qp, xc, UnitQuaternion(q..., false)), vector(qc))
    end

    # # contact origin points
    # cop = contact_point_origin(xp, qp, collision.contact_origin_parent) 
    # coc = contact_point_origin(xc, qc, collision.contact_origin_child)

    # # direction of minimum distance (child to parent)
    # d = cop - coc 
    # dir = normalize(d)

    # if relative == :parent 
    #     # cop - collision.contact_radius_parent * dir
    #     if jacobian == :parent 
    #         ∂c∂q = ∂contact_point_origin∂q(xp, qp, collision.contact_origin_parent)
    #         Q = ∂c∂q 
    #         Q -= collision.contact_radius_parent * ∂normalize∂x(d) * ∂c∂q
    #         return Q
    #     elseif jacobian == :child 
    #         Q = -1.0 * collision.contact_radius_parent * ∂normalize∂x(d) * -1.0 * ∂contact_point_origin∂q(xc, qc, collision.contact_origin_child)
    #         return Q
    #     end
    # elseif relative == :child 
    #     # coc + collision.contact_radius_child * dir
    #     if jacobian == :parent 
    #         Q = 1.0 * collision.contact_radius_child * ∂normalize∂x(d) * ∂contact_point_origin∂q(xp, qp, collision.contact_origin_parent)
    #         return Q
    #     elseif jacobian == :child 
    #         ∂c∂q = ∂contact_point_origin∂q(xc, qc, collision.contact_origin_child)
    #         Q = ∂c∂q 
    #         Q += collision.contact_radius_child * ∂normalize∂x(d) * -1.0 * ∂c∂q
    #         return Q
    #     end
    # end
end

# # normal projection (from child to parent)
# function contact_normal(collision::SphereSphereCollision, xp, qp, xc, qc)
#     # # contact origin points
#     cop = contact_point_origin(xp, qp, collision.contact_origin_parent) 
#     coc = contact_point_origin(xc, qc, collision.contact_origin_child)

#     # # unnormalized direction 
#     # dir0 = cop - coc

#     # # contact origin points
#     # cop = contact_point(:parent, collision, xp, qp, xc, qc) 
#     # coc = contact_point(:child,  collision, xp, qp, xc, qc)
 
#     # unnormalized direction 
#     dir = cop - coc

#     # distance 
#     dis = distance(collision, xp, qp, xc, qc) 

#     # normalized direction
#     return normalize(dir)'
# end

