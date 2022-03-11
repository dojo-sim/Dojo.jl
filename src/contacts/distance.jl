# function distance(collision::Collision, xp, qp, xc, qc)
#     collision.contact_normal * (xp + vector_rotate(collision.contact_origin, qp)) - collision.contact_radius
# end

# function ∂distance∂xp(collision::Collision, xp, qp, xc, qc)
#     collision.contact_normal
# end

# function ∂distance∂qp(collision::Collision, xp, qp, xc, qc)
#     collision.contact_normal * ∂vector_rotate∂q(collision.contact_origin, qp)
# end

