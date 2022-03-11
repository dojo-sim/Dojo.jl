"""
    ImpactContact{T,N} <: Contact{T,N}

    contact object for impact (i.e., no friction)

    collision: Collision
"""
mutable struct ImpactContact{T,N} <: Contact{T,N}
    collision::Collision{T,0,3,0}

    function ImpactContact(body::Body{T}, normal::AbstractVector{T}; 
        contact_origin=szeros(T, 3), 
        contact_radius=0.0) where T
        # projector
        V1, V2, V3 = orthogonal_columns(normal) #
        A = [V1 V2 V3]
        Ainv = inv(A)
        contact_normal = Ainv[3, SA[1; 2; 3]]'
        # collision 
        collision = SphereFloorCollision(szeros(T, 0, 3), contact_normal, SVector{3}(contact_origin), contact_radius)
        new{Float64,2}(collision)
    end
end

function constraint(mechanism, contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:ImpactContact{T,N}}
    # contact model
    model = contact.model

    # parent
    pbody = get_body(mechanism, contact.parent_id)
    xp, qp = next_configuration(pbody.state, mechanism.timestep)

    # child
    cbody = get_body(mechanism, contact.child_id)
    xc, qc = next_configuration(cbody.state, mechanism.timestep)

    # constraint 
    SVector{1,T}(distance(model.collision, xp, qp, xc, qc) - contact.impulses_dual[2][1])
end

function constraint_jacobian(contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:ImpactContact{T,N},N½}
    γ = contact.impulses[2] + REG * neutral_vector(contact.model)
    s = contact.impulses_dual[2] + REG * neutral_vector(contact.model)
    ∇s = hcat(γ, -Diagonal(sones(N½)))
    ∇γ = hcat(s, -Diagonal(szeros(N½)))
    return [∇s ∇γ]
end

function constraint_jacobian_configuration(model::ImpactContact, 
    xp::AbstractVector, vp::AbstractVector, qp::UnitQuaternion, ϕp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::UnitQuaternion, ϕc::AbstractVector, 
    timestep)

    X = ∂distance∂xp(model.collision, xp, qp, xc, qc)
    Q = ∂distance∂qp(model.collision, xp, qp, xc, qc)

    return [X Q]
end

function constraint_jacobian_velocity(model::ImpactContact, 
    xp::AbstractVector, vp::AbstractVector, qp::UnitQuaternion, ϕp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::UnitQuaternion, ϕc::AbstractVector, 
    timestep)

    # recover current orientation 
    qp_prev = next_orientation(qp, -ϕp, timestep)

    # Jacobian
    V = ∂distance∂xp(model.collision, xp, qp, xc, qc) * timestep
    Ω = ∂distance∂qp(model.collision, xp, qp, xc, qc) * rotational_integrator_jacobian_velocity(qp_prev, ϕp, timestep)

    return [V Ω]
end

function force_mapping(model::ImpactContact, xp::AbstractVector, qp::UnitQuaternion, xc::AbstractVector, qc::UnitQuaternion)
    X = model.collision.contact_normal
    return X
end

