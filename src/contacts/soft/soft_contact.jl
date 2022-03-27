
"""
    SoftContact{T,N} <: Contact{T,N}

    contact object for soft contact

    collision: Collision
"""
mutable struct SoftContact{T,N} <: Contact{T,N}
    friction_parameterization::SMatrix{2,2,T,4}
    collision::SoftCollision{T,2,3,6}
    contact_origin::SVector{3,T} # origin of the frame in which the nerf is defined
end

function SoftContact(body::Body{T}, normal::AbstractVector{T}, collider::Collider;
    contact_origin=-collider.center_of_mass,
    contact_radius=0.0) where T

    # contact directions
    V1, V2, V3 = orthogonal_columns(normal) #
    A = [V1 V2 V3]
    Ainv = inv(A)
    contact_normal = Ainv[3, SA[1; 2; 3]]'
    contact_tangent = Ainv[SA[1; 2], SA[1; 2; 3]]

    # friction parametrization
    parameterization = SA{T}[
         1.0  0.0
         0.0  1.0
    ]

    # collision
    collision = SoftHalfSpaceCollision(collider, contact_tangent, contact_normal, SVector{3}(contact_origin))

    SoftContact{Float64,6}(parameterization, collision, contact_origin)
end

function constraint(mechanism, contact::SoftContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:SoftContact{T,N}}
    timestep = mechanism.timestep
    pbody = get_body(mechanism, contact.parent_id)
    cbody = get_body(mechanism, contact.child_id)
    xp, vp, qp, ϕp = next_configuration_velocity(pbody.state, timestep)
    xc, vc, qc, ϕc = next_configuration_velocity(cbody.state, timestep)
    impulse = constraint(contact.model, xp, vp, qp, ϕp, xc, vc, qc, ϕc, timestep)

    VERBOSE && println(
        # "x2p:", scn.(x2p),
        # "    ψ:", scn(ψ),
        "    λ:", scn.(contact.impulses[2]),
        "    imp:", scn.(impulse),
        # "    err:", scn(impulse[1] - contact.impulses[2][1])
        )
    return impulse - contact.impulses[2]
end

function constraint_jacobian(contact::SoftContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:SoftContact{T,N}}
    ∇γ = -Diagonal(sones(N))
    return ∇γ
end

function constraint_jacobian_configuration(relative::Symbol, model::SoftContact{T,N},
    xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ϕp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ϕc::AbstractVector,
    timestep) where {T,N}
    return szeros(T,N,6)
end

function constraint_jacobian_velocity(relative::Symbol, model::SoftContact{T,N},
    xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ϕp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ϕc::AbstractVector,
    timestep) where {T,N}

    collider = model.collision.collider
    opts = collider.options
    smoothing = opts.coulomb_smoothing
    regularizer = opts.coulomb_regularizer
    ψ, barycenter, contact_normal, vc_normal, vc_tangential, F_impact, F_friction, impulse = constraint(model,
        xp, vp, qp, ϕp,
        xc, vc, qc, ϕc,
        timestep, intermediate=true)
    # recover current orientation
    q2p = next_orientation(qp, -ϕp, timestep)
    ∂vc∂v = Diagonal(sones(3))
    ∂vc∂ϕ = rotation_matrix(q2p) * skew(collider.center_of_mass - barycenter)
    if relative == :parent
        # impact
        ∂F_impact∂vc = -ψ * opts.impact_damper * contact_normal * contact_normal'
        FV_impact = ∂F_impact∂vc * ∂vc∂v
        FΩ_impact = ∂F_impact∂vc * ∂vc∂ϕ
        # friction
        ∂F_friction∂vc =  -opts.sliding_drag * vc_tangential * (FV_impact * F_impact/(norm(F_impact) + 1e-20))'
        ∂F_friction∂vc += -opts.sliding_drag * norm(F_impact) * (Diagonal(sones(3)) - contact_normal * contact_normal')
        ∂F_friction∂vc +=  opts.sliding_friction *
            coulomb_direction(vc_tangential, smoothing, regularizer) *
            (FV_impact * F_impact/(norm(F_impact) + 1e-20))'
        ∂F_friction∂vc +=  opts.sliding_friction * norm(F_impact) *
            ∂coulomb_direction∂v(vc_tangential, smoothing, regularizer) *
            (Diagonal(sones(3)) - contact_normal * contact_normal')
        FV_friction = ∂F_friction∂vc * ∂vc∂v
        FΩ_friction = ∂F_friction∂vc * ∂vc∂ϕ

        τV_contact = (-opts.rolling_drag * ϕp + opts.rolling_friction * coulomb_direction(ϕp, smoothing, regularizer)) *
            (FV_impact * F_impact/(norm(F_impact) + 1e-20))'
        τΩ_contact = τV_contact * ∂vc∂ϕ
        τΩ_contact += norm(F_impact) * (-opts.rolling_drag * Diagonal(sones(3)) +
            opts.rolling_friction * ∂coulomb_direction∂v(ϕp, smoothing, regularizer))
        # contact
        V = [FV_impact + FV_friction; τV_contact]
        Ω = [FΩ_impact + FΩ_friction; τΩ_contact]
    elseif relative == :child
        V = szeros(T,N,3)
        Ω = szeros(T,N,3)
    end
    return timestep * [V Ω]
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
