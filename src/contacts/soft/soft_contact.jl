
"""
    SoftContact{T,N} <: Contact{T,N}

    contact object for soft contact

    collision: Collision
"""
mutable struct SoftContact{T,N,C} <: Contact{T,N}
    collision::C
    collider_origin::SVector{3,T} # origin of the frame in which the nerf is defined
    ψ::T # amount of violation
    barycenter::SVector{3,T} # position of the overlap region's barycenter in the parent frame (=collider frame)
    normal::SVector{3,T} # contact normal expressed in the world frame
end

function SoftContact(body::Body{T}, normal::AbstractVector{T}, collider::Collider;
        collider_origin=-collider.center_of_mass, sphere_origin=szeros(T,3), radius::T=0.0,
        collision_type::Symbol=:soft_halfspace) where T

    # contact directions
    V1, V2, V3 = orthogonal_columns(normal)
    A = [V1 V2 V3]
    Ainv = inv(A)
    contact_normal = Ainv[3, SA[1; 2; 3]]'
    contact_tangent = Ainv[SA[1; 2], SA[1; 2; 3]]

    # collision
    if collision_type == :soft_halfspace
        collision = SoftHalfSpaceCollision(collider, contact_tangent, contact_normal, SVector{3}(collider_origin))
    elseif collision_type == :soft_sphere
        collision = SoftSphereCollision(collider, contact_tangent, SVector{3}(collider_origin), sphere_origin, radius)
    else
        error("Unknown collision_type")
    end
    SoftContact{Float64,6,typeof(collision)}(collision, collider_origin, 0.0, szeros(3), szeros(3))
end

function constraint(mechanism, contact::SoftContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:SoftContact{T,N}}
    timestep = mechanism.timestep
    pbody = get_body(mechanism, contact.parent_id)
    cbody = get_body(mechanism, contact.child_id)
    xp, vp, qp, ϕp = next_configuration_velocity(pbody.state, timestep)
    xc, vc, qc, ϕc = next_configuration_velocity(cbody.state, timestep)
    impulse = constraint(contact.model, xp, vp, qp, ϕp, xc, vc, qc, ϕc, timestep)
    return impulse - contact.impulses[2]
end

function constraint_jacobian(contact::SoftContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:SoftContact{T,N}}
    ∇γ = -Diagonal(sones(N))
    return ∇γ
end

function constraint(model::SoftContact{T},
        xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ϕp::AbstractVector,
        xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ϕc::AbstractVector,
        timestep; intermediate::Bool=false) where T

    collision = model.collision
    collider = collision.collider
    opts = collider.options
    smoothing = opts.coulomb_smoothing
    regularizer = opts.coulomb_regularizer
    x2p = next_position(xp, -vp, timestep)
    q2p = next_orientation(qp, -ϕp, timestep)
    x2c = next_position(xc, -vc, timestep)
    q2c = next_orientation(qc, -ϕc, timestep)

    # collision
    # barycenter in soft_body frame
    # normal in world frame
    ψ, barycenter, normal = model.ψ, model.barycenter, model.normal
    # velocities
    vp_ = vp + vector_rotate(skew(collider.center_of_mass - barycenter) * ϕp, q2p)
    xp_ = x2p + Dojo.vector_rotate(barycenter - collider.center_of_mass, q2p)
    vc_ = vc + skew(x2c - xp_) * vector_rotate(ϕc, q2c)
    v = vp_ - vc_
    v_normal = normal * normal' * v
    v_tangential = v - v_normal
    # impact
    F_impact = ψ * opts.impact_spring * normal
    F_impact -= ψ * opts.impact_damper * v_normal
    # friction
    F_friction = -opts.sliding_drag * norm(F_impact) * v_tangential
    F_friction += opts.sliding_friction * norm(F_impact) * coulomb_direction(v_tangential, smoothing, regularizer)
    F_contact = F_impact + F_friction

    ϕ = vector_rotate(ϕp, q2p) - vector_rotate(ϕc, q2c)
    τ_contact = -opts.rolling_drag * norm(F_impact) * ϕ
    τ_contact += opts.rolling_friction * norm(F_impact) * coulomb_direction(ϕ, smoothing, regularizer)

    # constraint
    impulse = timestep * [F_contact; τ_contact]
    !intermediate && return impulse
    return ψ, barycenter, normal, v_normal, v_tangential, F_impact, F_friction, impulse
end

function constraint_jacobian_velocity(relative::Symbol, model::SoftContact{T,N},
    xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ϕp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ϕc::AbstractVector,
    timestep) where {T,N}

    impulse = constraint(model, xp, vp, qp, ϕp, xc, vc, qc, ϕc, timestep)
    if relative == :parent
        ∂impulse∂vϕ = FiniteDiff.finite_difference_jacobian(
            vϕ -> constraint(model, xp, vϕ[SUnitRange(1,3)], qp, vϕ[SUnitRange(4,6)], xc, vc, qc, ϕc, timestep),
            [vp; ϕp])
    elseif relative == :child
        ∂impulse∂vϕ = FiniteDiff.finite_difference_jacobian(
            vϕ -> constraint(model, xp, vp, qp, ϕp, xc, vϕ[SUnitRange(1,3)], qc, vϕ[SUnitRange(4,6)], timestep),
            [vc; ϕc])
    end
    return ∂impulse∂vϕ
end

function constraint_jacobian_configuration(relative::Symbol, model::SoftContact{T,N},
    xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ϕp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ϕc::AbstractVector,
    timestep) where {T,N}
    return szeros(T,N,6)
end

function force_mapping(relative::Symbol, model::SoftContact,
    xp::AbstractVector, qp::Quaternion,
    xc::AbstractVector, qc::Quaternion)

    X = Diagonal(sones(3)) # TODO hard-coded
    if relative == :parent
        return X
    elseif relative == :child
        return -1.0 * X
    end
end

function ∂force_mapping_jvp∂x(relative::Symbol, jacobian::Symbol,
    model::SoftContact{T},
    xp::AbstractVector, qp::Quaternion,
    xc::AbstractVector, qc::Quaternion,
    λ::AbstractVector) where T

    X = szeros(T,3,3)
    if relative == :parent
        return X
    elseif relative == :child
        return -1.0 * X
    end
end

function ∂force_mapping_jvp∂q(relative::Symbol, jacobian::Symbol,
    model::SoftContact{T},
    xp::AbstractVector, qp::Quaternion,
    xc::AbstractVector, qc::Quaternion,
    λ::AbstractVector) where T

    X = szeros(T,3,4)

    if relative == :parent
        return X
    elseif relative == :child
        return -1.0 * X
    end
end
