# contact point velocity
function contact_point_velocity(model::Contact, x, q, v, ω, c)
    return v + skew(vector_rotate(ω, q)) * (c - x)
end

function ∂contact_point_velocity∂x(model::Contact, x, q, v, ω, c)
    return -skew(vector_rotate(ω, q))
end

function ∂contact_point_velocity∂q(model::Contact, x, q, v, ω, c)
    return -skew(c - x) * ∂vector_rotate∂q(ω, q)
end

function ∂contact_point_velocity∂v(model::Contact, x, q, v, ω, c)
    return 1.0 * sI(3)
end

function ∂contact_point_velocity∂ω(model::Contact, x, q, v, ω, c)
    return -skew(c - x) * rotation_matrix(q)
end

function ∂contact_point_velocity∂c(model::Contact, x, q, v, ω, c)
    return skew(vector_rotate(ω, q))
end

# relative tangential velocity
function relative_tangential_velocity(model::Contact, xp, qp, vp, ωp, xc, qc, vc, ωc)
    # contact points
    cp = contact_point(:parent, model.collision, xp, qp, xc, qc)
    cc = contact_point(:child, model.collision, xp, qp, xc, qc)

    # contact point velocities
    vp = contact_point_velocity(model, xp, qp, vp, ωp, cp)
    vc = contact_point_velocity(model, xc, qc, vc, ωc, cc)

    # tangent projection
    return contact_tangent(model.collision, xp, qp, xc, qc) * (vp - vc)
end

function ∂relative_tangential_velocity∂x(jacobian::Symbol, model::Contact, xp, qp, vp, ωp, xc, qc, vc, ωc)
    # contact points
    cp = contact_point(:parent, model.collision, xp, qp, xc, qc)
    cc = contact_point(:child, model.collision, xp, qp, xc, qc)

    # contact point velocities
    vp = contact_point_velocity(model, xp, qp, vp, ωp, cp)
    vc = contact_point_velocity(model, xc, qc, vc, ωc, cc)

    Δv = vp - vc

    # Jacobians
    if jacobian == :parent
        X = contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂x(model, xp, qp, vp, ωp, cp)
        X += contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xp, qp, vp, ωp, cp) * ∂contact_point∂x(:parent, jacobian, model.collision, xp, qp, xc, qc)
        X -= contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xc, qc, vc, ωc, cc) * ∂contact_point∂x(:child, jacobian, model.collision, xp, qp, xc, qc)
    elseif jacobian == :child
        X = contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xp, qp, vp, ωp, cp) * ∂contact_point∂x(:parent, jacobian, model.collision, xp, qp, xc, qc)
        X -= contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂x(model, xc, qc, vc, ωc, cc)
        X -= contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xc, qc, vc, ωc, cc) * ∂contact_point∂x(:child, jacobian, model.collision, xp, qp, xc, qc)
    end

    X += [
        Δv' * ∂contact_tangent_one_transpose∂x(jacobian, model.collision, xp, qp, xc, qc);
        Δv' * ∂contact_tangent_two_transpose∂x(jacobian, model.collision, xp, qp, xc, qc);
    ]

    return X
end

function ∂relative_tangential_velocity∂q(jacobian::Symbol, model::Contact, xp, qp, vp, ωp, xc, qc, vc, ωc)
    # contact points
    cp = contact_point(:parent, model.collision, xp, qp, xc, qc)
    cc = contact_point(:child, model.collision, xp, qp, xc, qc)

    # contact point velocities
    vp = contact_point_velocity(model, xp, qp, vp, ωp, cp)
    vc = contact_point_velocity(model, xc, qc, vc, ωc, cc)

    Δv = vp - vc

    # Jacobians
    if jacobian == :parent
        X = contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂q(model, xp, qp, vp, ωp, cp)
        X += contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xp, qp, vp, ωp, cp) * ∂contact_point∂q(:parent, jacobian, model.collision, xp, qp, xc, qc)
        X -= contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xc, qc, vc, ωc, cc) * ∂contact_point∂q(:child, jacobian, model.collision, xp, qp, xc, qc)
    elseif jacobian == :child
        X = contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xp, qp, vp, ωp, cp) * ∂contact_point∂q(:parent, jacobian, model.collision, xp, qp, xc, qc)
        X -= contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂q(model, xc, qc, vc, ωc, cc)
        X -= contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂c(model, xc, qc, vc, ωc, cc) * ∂contact_point∂q(:child, jacobian, model.collision, xp, qp, xc, qc)
    end

    # X += ∂contact_tangent_jvp∂q(jacobian, model.collision, xp, qp, xc, qc, Δv)
    X += [
        Δv' * ∂contact_tangent_one_transpose∂q(jacobian, model.collision, xp, qp, xc, qc);
        Δv' * ∂contact_tangent_two_transpose∂q(jacobian, model.collision, xp, qp, xc, qc);
    ]

    return X
end

function ∂relative_tangential_velocity∂v(jacobian::Symbol, model::Contact, xp, qp, vp, ωp, xc, qc, vc, ωc)
    # contact points
    cp = contact_point(:parent, model.collision, xp, qp, xc, qc)
    cc = contact_point(:child, model.collision, xp, qp, xc, qc)

    # contact point velocities
    vp = contact_point_velocity(model, xp, qp, vp, ωp, cp)
    vc = contact_point_velocity(model, xc, qc, vc, ωc, cc)

    # Δv = vp - vc

    # Jacobians
    if jacobian == :parent
        X = contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂v(model, xp, qp, vp, ωp, cp)
    elseif jacobian == :child
        X = -1.0 * contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂v(model, xc, qc, vc, ωc, cc)
    end

    return X
end

function ∂relative_tangential_velocity∂ω(jacobian::Symbol, model::Contact, xp, qp, vp, ωp, xc, qc, vc, ωc)
    # contact points
    cp = contact_point(:parent, model.collision, xp, qp, xc, qc)
    cc = contact_point(:child, model.collision, xp, qp, xc, qc)

    # contact point velocities
    vp = contact_point_velocity(model, xp, qp, vp, ωp, cp)
    vc = contact_point_velocity(model, xc, qc, vc, ωc, cc)

    # Δv = vp - vc

    if jacobian == :parent
        X = contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂ω(model, xp, qp, vp, ωp, cp)
    elseif jacobian == :child
        X = -1.0 * contact_tangent(model.collision, xp, qp, xc, qc) * ∂contact_point_velocity∂ω(model, xc, qc, vc, ωc, cc)
    end

    return X
end
