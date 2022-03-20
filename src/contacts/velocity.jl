# contact point velocity
function contact_point_velocity(model::Contact, x, q, v, ϕ, c)
    return v + skew(vector_rotate(ϕ, q)) * (c - x)
end

function ∂contact_point_velocity∂x(model::Contact, x, q, v, ϕ, c)
    return -skew(vector_rotate(ϕ, q))
end

function ∂contact_point_velocity∂q(model::Contact, x, q, v, ϕ, c)
    return -skew(c - x) * ∂vector_rotate∂q(ϕ, q)
end

function ∂contact_point_velocity∂v(model::Contact, x, q, v, ϕ, c)
    return 1.0 * I(3)
end

function ∂contact_point_velocity∂ϕ(model::Contact, x, q, v, ϕ, c)
    return -skew(c - x) * ∂vector_rotate∂p(ϕ, q)
end

function ∂contact_point_velocity∂c(model::Contact, x, q, v, ϕ, c)
    return skew(vector_rotate(ϕ, q))
end

# relative tangential velocity
function relative_tangential_velocity(model::Contact, xp, qp, vp, ϕp, xc, qc, vc, ϕc)
    # contact points
    cp = contact_point(:parent, model.collision, xp, qp, xc, qc)
    cc = contact_point(:child, model.collision, xp, qp, xc, qc)

    # contact point velocities
    vp = contact_point_velocity(model, xp, qp, vp, ϕp, cp)
    vc = contact_point_velocity(model, xc, qc, vc, ϕc, cc)

    # tangent projection
    return contact_tangent(model.collision, xp, qp, xc, qc) * (vp - vc)
end

function ∂relative_tangential_velocity∂x(jacobian::Symbol, model::Contact, xp, qp, vp, ϕp, xc, qc, vc, ϕc)
    # contact points
    cp = contact_point(:parent, model.collision, xp, qp, xc, qc)
    cc = contact_point(:child, model.collision, xp, qp, xc, qc)

    # contact point velocities
    vp = contact_point_velocity(model, xp, qp, vp, ϕp, cp)
    vc = contact_point_velocity(model, xc, qc, vc, ϕc, cc)

    Δv = vp - vc

    # Jacobians 
    if jacobian == :parent
        X = contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂x(model, xp, qp, vp, ϕp, cp) 
        X += contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xp, qp, vp, ϕp, cp) * ∂contact_point∂x(:parent, jacobian, model.collision, xp, qp, xc, qc)
        X -= contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xc, qc, vc, ϕc, cc) * ∂contact_point∂x(:child, jacobian, model.collision, xp, qp, xc, qc)
    elseif jacobian == :child 
        X = contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xp, qp, vp, ϕp, cp) * ∂contact_point∂x(:parent, jacobian, model.collision, xp, qp, xc, qc)
        X -= contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂x(model, xc, qc, vc, ϕc, cc) 
        X -= contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xc, qc, vc, ϕc, cc) * ∂contact_point∂x(:child, jacobian, model.collision, xp, qp, xc, qc)
    end

    X += [
        Δv' * ∂contact_tangent_one_tangent∂x(jacobian, model.collision, xp, qp, xc, qc);
        Δv' * ∂contact_tangent_two_tangent∂x(jacobian, model.collision, xp, qp, xc, qc);
    ]

    return X
end

function ∂relative_tangential_velocity∂q(jacobian::Symbol, model::Contact, xp, qp, vp, ϕp, xc, qc, vc, ϕc)
    # contact points
    cp = contact_point(:parent, model.collision, xp, qp, xc, qc)
    cc = contact_point(:child, model.collision, xp, qp, xc, qc)

    # contact point velocities
    vp = contact_point_velocity(model, xp, qp, vp, ϕp, cp)
    vc = contact_point_velocity(model, xc, qc, vc, ϕc, cc)

    Δv = vp - vc

    # Jacobians 
    if jacobian == :parent
        X = contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂q(model, xp, qp, vp, ϕp, cp) 
        X += contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xp, qp, vp, ϕp, cp) * ∂contact_point∂q(:parent, jacobian, model.collision, xp, qp, xc, qc)
        X -= contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xc, qc, vc, ϕc, cc) * ∂contact_point∂q(:child, jacobian, model.collision, xp, qp, xc, qc)
    elseif jacobian == :child 
        X = contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xp, qp, vp, ϕp, cp) * ∂contact_point∂q(:parent, jacobian, model.collision, xp, qp, xc, qc)
        X -= contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂q(model, xc, qc, vc, ϕc, cc) 
        X -= contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xc, qc, vc, ϕc, cc) * ∂contact_point∂q(:child, jacobian, model.collision, xp, qp, xc, qc)
    end

    # X += ∂contact_tangent_jvp∂q(jacobian, model.collision, xp, qp, xc, qc, Δv)
    X += [
        Δv' * ∂contact_tangent_one_transpose∂q(jacobian, model.collision, xp, qp, xc, qc);
        Δv' * ∂contact_tangent_two_transpose∂q(jacobian, model.collision, xp, qp, xc, qc);
    ]

    return X
end

function ∂relative_tangential_velocity∂v(jacobian::Symbol, model::Contact, xp, qp, vp, ϕp, xc, qc, vc, ϕc)
    # contact points
    cp = contact_point(:parent, model.collision, xp, qp, xc, qc)
    cc = contact_point(:child, model.collision, xp, qp, xc, qc)

    # contact point velocities
    vp = contact_point_velocity(model, xp, qp, vp, ϕp, cp)
    vc = contact_point_velocity(model, xc, qc, vc, ϕc, cc)

    # Δv = vp - vc

    # Jacobians 
    if jacobian == :parent
        X = contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂v(model, xp, qp, vp, ϕp, cp) 
    elseif jacobian == :child 
        X = -1.0 * contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂v(model, xc, qc, vc, ϕc, cc) 
    end

    return X
end

function ∂relative_tangential_velocity∂ϕ(jacobian::Symbol, model::Contact, xp, qp, vp, ϕp, xc, qc, vc, ϕc)
    # contact points
    cp = contact_point(:parent, model.collision, xp, qp, xc, qc)
    cc = contact_point(:child, model.collision, xp, qp, xc, qc)

    # contact point velocities
    vp = contact_point_velocity(model, xp, qp, vp, ϕp, cp)
    vc = contact_point_velocity(model, xc, qc, vc, ϕc, cc)

    # Δv = vp - vc

    if jacobian == :parent
        X = contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂ϕ(model, xp, qp, vp, ϕp, cp) 
    elseif jacobian == :child 
        X = -1.0 * contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂ϕ(model, xc, qc, vc, ϕc, cc) 
    end

    return X
end



