"""
    SoftHalfSpaceCollision

    collision between a soft contact and a flat surface

    collider: soft object
    contact_tangent: mapping from world frame to surface tangent frame
    contact_normal: inverse/complement of contact_tangent
    contact_origin: position of contact on Body relative to center of mass
"""
mutable struct SoftHalfSpaceCollision{T,O,I,OI,N} <: SoftCollision{T,O,I,OI,N}
    collider::SoftCollider{T,N}
    contact_tangent::SMatrix{O,I,T,OI}
    contact_normal::Adjoint{T,SVector{I,T}}
    contact_origin::SVector{I,T}
end

function inside(collision::SoftHalfSpaceCollision{T}, p) where T
    origin = szeros(T,3)
    normal = collision.contact_normal
    c = normal * (p - origin)
    return c <= 0.0
end

# normal projection (from child to parent)
function contact_normal(collision::SoftHalfSpaceCollision, p, xc, qc)
    return collision.contact_normal[1,:]
end

#
#
# # distance
# function distance(collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     collision.contact_normal * (xp + vector_rotate(collision.contact_origin, qp)) - collision.contact_radius
# end
#
# function ∂distance∂x(gradient::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     if gradient == :parent
#         return collision.contact_normal
#     elseif gradient == :child
#         return szeros(eltype(xc), 1, 3)
#     end
# end
#
# function ∂distance∂q(gradient::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     if gradient == :parent
#         return collision.contact_normal * ∂vector_rotate∂q(collision.contact_origin, qp)
#     elseif gradient == :child
#         return szeros(eltype(qc), 1, 4)
#     end
# end
#
# # contact point in world frame
# function contact_point(relative::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     if relative == :parent
#         return xp + vector_rotate(collision.contact_origin, qp) - collision.contact_normal' * collision.contact_radius
#     elseif relative == :child
#         projector = collision.contact_tangent' * collision.contact_tangent
#         return projector * (xp + vector_rotate(collision.contact_origin, qp))
#     end
# end
#
# function ∂contact_point∂x(relative::Symbol, jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     if relative == :parent
#         if jacobian == :parent
#             return 1.0 * I(3)
#         elseif jacobian == :child
#             return szeros(eltype(xp), 3, 3)
#         end
#     elseif relative == :child
#         if jacobian == :parent
#             projector = collision.contact_tangent' * collision.contact_tangent
#             return projector
#         elseif jacobian == :child
#             return szeros(eltype(xp), 3, 3)
#         end
#     end
# end
#
# function ∂contact_point∂q(relative::Symbol, jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     if relative == :parent
#         if jacobian == :parent
#             return ∂vector_rotate∂q(collision.contact_origin, qp)
#         end
#     elseif relative == :child
#         if jacobian == :parent
#             projector = collision.contact_tangent' * collision.contact_tangent
#             return projector * ∂vector_rotate∂q(collision.contact_origin, qp)
#         elseif jacobian == :child
#             return szeros(eltype(qp), 3, 4)
#         end
#     end
# end
#
# # normal projection (from child to parent)
# function contact_normal(collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return collision.contact_normal
# end
#
# function ∂contact_normal_transpose∂x(jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return szeros(eltype(collision.contact_normal), 3, 3)
# end
#
# function ∂contact_normal_transpose∂q(jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return szeros(eltype(collision.contact_normal), 3, 4)
# end
#
# # tangent projection
# function contact_tangent(collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return collision.contact_tangent # {2,4} x 3
# end
#
# function ∂contact_tangent_one_transpose∂x(jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return szeros(eltype(collision.contact_tangent), 3, 3)
# end
#
# function ∂contact_tangent_two_transpose∂x(jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return szeros(eltype(collision.contact_tangent), 3, 3)
# end
#
# function ∂contact_tangent_one_transpose∂q(jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return szeros(eltype(collision.contact_tangent), 3, 4)
# end
#
# function ∂contact_tangent_two_transpose∂q(jacobian::Symbol, collision::SoftHalfSpaceCollision, xp, qp, xc, qc)
#     return szeros(eltype(collision.contact_tangent), 3, 4)
# end
