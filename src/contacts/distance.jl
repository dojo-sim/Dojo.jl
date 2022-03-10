function distance(collision::Collision, xp, qp, xc, qc)
    collision.surface_normal_projector * (xp + vector_rotate(collision.contact_point, qp)) - collision.contact_radius
end

function ∂distance∂xp(collision::Collision, xp, qp, xc, qc)
    collision.surface_normal_projector
end

function ∂distance∂qp(collision::Collision, xp, qp, xc, qc)
    collision.surface_normal_projector * ∂vector_rotate∂q(collision.contact_point, qp)
end

