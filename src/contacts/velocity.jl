# contact point velocity
function contact_point_velocity(model::Contact, x, q, v, ϕ)
    return v + skew(vector_rotate(ϕ, q)) * (vector_rotate(model.collision.contact_point, q) - model.collision.surface_normal_projector' * model.collision.contact_radius)
end

function ∂contact_point_velocity∂x(model::Contact, x, q, v, ϕ)
    return szeros(3, 3)
end

function ∂contact_point_velocity∂q(model::Contact, x, q, v, ϕ)
    ∂v∂q = skew(vector_rotate(ϕ, q)) * ∂vector_rotate∂q(model.collision.contact_point, q)
    ∂v∂q += skew(model.collision.surface_normal_projector' * model.collision.contact_radius - vector_rotate(model.collision.contact_point, q)) * ∂vector_rotate∂q(ϕ, q)
    return ∂v∂q
end

function ∂contact_point_velocity∂v(model::Contact, x, q, v, ϕ)
    return 1.0 * I(3)
end

function ∂contact_point_velocity∂ϕ(model::Contact, x, q, v, ϕ)
    ∂v∂ϕ = skew(model.collision.surface_normal_projector' * model.collision.contact_radius - vector_rotate(model.collision.contact_point, q)) * ∂vector_rotate∂p(ϕ, q)
    return ∂v∂ϕ
end

# relative tangential velocity
function relative_tangential_velocity(model::Contact, xp, qp, vp, ϕp, xc, qc, vc, ϕc)
    vp = contact_point_velocity(model, xp, qp, vp, ϕp)
    vc = contact_point_velocity(model, xc, qc, vc, ϕc)
    return model.collision.surface_projector * (vp - vc)
end

function ∂relative_tangential_velocity∂x(jacobian::Symbol, model::Contact, xp, qp, vp, ϕp, xc, qc, vc, ϕc)
    if jacobian == :parent 
        return model.collision.surface_projector * ∂contact_point_velocity∂x(model, xp, qp, vp, ϕp)
    elseif jacobian == :child 
        return model.collision.surface_projector * -1.0 * ∂contact_point_velocity∂x(model, xc, qc, vc, ϕc)
    end
end

function ∂relative_tangential_velocity∂q(jacobian::Symbol, model::Contact, xp, qp, vp, ϕp, xc, qc, vc, ϕc)
    if jacobian == :parent 
        return model.collision.surface_projector * ∂contact_point_velocity∂q(model, xp, qp, vp, ϕp)
    elseif jacobian == :child 
        return model.collision.surface_projector * -1.0 * ∂contact_point_velocity∂q(model, xc, qc, vc, ϕc)
    end
end

function ∂relative_tangential_velocity∂v(jacobian::Symbol, model::Contact, xp, qp, vp, ϕp, xc, qc, vc, ϕc)
    if jacobian == :parent 
        return model.collision.surface_projector * ∂contact_point_velocity∂v(model, xp, qp, vp, ϕp)
    elseif jacobian == :child 
        return model.collision.surface_projector * -1.0 * ∂contact_point_velocity∂v(model, xc, qc, vc, ϕc)
    end
end

function ∂relative_tangential_velocity∂ϕ(jacobian::Symbol, model::Contact, xp, qp, vp, ϕp, xc, qc, vc, ϕc)
    if jacobian == :parent 
        return model.collision.surface_projector * ∂contact_point_velocity∂ϕ(model, xp, qp, vp, ϕp)
    elseif jacobian == :child 
        return model.collision.surface_projector * -1.0 * ∂contact_point_velocity∂ϕ(model, xc, qc, vc, ϕc)
    end
end


