function distance(model::Contact, xp, qp, xc, qc)
    model.surface_normal_projector * (xp + vector_rotate(model.contact_point, qp) - model.offset)
end

function ∂distance∂xp(model::Contact, xp, qp, xc, qc)
    model.surface_normal_projector
end

function ∂distance∂qp(model::Contact, xp, qp, xc, qc)
    model.surface_normal_projector * ∂vector_rotate∂q(model.contact_point, qp)
end

