
"""
    SoftContact{T,N} <: Contact{T,N}

    contact object for soft contact

    collision: Collision
"""
mutable struct SoftContact{T,N,C} <: Contact{T,N}
    collision::C
    ψ::T # amount of violation
    barycenter::SVector{3,T} # position of the overlap region's barycenter in the parent frame (=collider frame)
    normal::SVector{3,T} # contact normal expressed in the world frame
end

function SoftContact(body::Body{T}, normal::AbstractVector{T}, collider::Collider;
        parent_origin=-collider.center_of_mass, child_origin=szeros(T,3), radius::T=0.0,
        child_collider=nothing, child_collider_origin=nothing,
        collision_type::Symbol=:soft_halfspace) where T

    if child_collider != nothing && child_collider_origin == nothing
        child_collider_origin = -child_collider.center_of_mass
    end
    # contact directions
    V1, V2, V3 = orthogonal_columns(normal)
    A = [V1 V2 V3]
    Ainv = inv(A)
    contact_normal = Ainv[3, SA[1; 2; 3]]'
    contact_tangent = Ainv[SA[1; 2], SA[1; 2; 3]]

    # collision
    if collision_type == :soft_halfspace
        collision = SoftHalfSpaceCollision(collider, contact_normal, SVector{3}(parent_origin), contact_tangent)
    elseif collision_type == :soft_sphere
        collision = SoftSphereCollision(collider, SVector{3}(parent_origin), child_origin, radius, contact_tangent)
    elseif collision_type == :soft_soft
        collision = SoftSoftCollision(collider, child_collider, SVector{3}(parent_origin), SVector{3}(child_collider_origin), contact_tangent)
    else
        error("Unknown collision_type")
    end
    # SoftContact{Float64,6,typeof(collision)}(collision, collider_origin, 0.0, szeros(3), szeros(3))
    SoftContact{Float64,6,typeof(collision)}(collision, 0.0, szeros(3), szeros(3))
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
        timestep; recompute::Bool=false) where T

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
    if recompute
        model.ψ, model.barycenter, model.normal = overlap(collision, x2p, q2p, x2c, q2c)
    end
    ψ, barycenter, normal = model.ψ, model.barycenter, model.normal
    # velocities
    constraint(collider, ψ, barycenter, normal, x2p, vp, q2p, ϕp, x2c, vc, q2c, ϕc, timestep)
end


function constraint(collider, ψ, barycenter, normal,
        x2p::AbstractVector, vp::AbstractVector, q2p::Quaternion, ϕp::AbstractVector,
        x2c::AbstractVector, vc::AbstractVector, q2c::Quaternion, ϕc::AbstractVector,
        timestep) where T

    opts = collider.options
    smoothing = opts.coulomb_smoothing
    regularizer = opts.coulomb_regularizer

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
    return impulse
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
    # @show ∂impulse∂vϕ[:,1]
    return ∂impulse∂vϕ
end

function constraint_jacobian_configuration(relative::Symbol, model::SoftContact{T,N},
    xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ϕp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ϕc::AbstractVector,
    timestep) where {T,N}
    @show "this might be wrong I am not sure wrt which configuration we are diffing"
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
