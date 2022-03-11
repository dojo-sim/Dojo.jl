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
        collision = SphereFlatCollision(szeros(T, 0, 3), contact_normal, SVector{3}(contact_origin), contact_radius)

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

function constraint_jacobian_configuration(relative::Symbol, model::ImpactContact, 
    xp::AbstractVector, vp::AbstractVector, qp::UnitQuaternion, ϕp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::UnitQuaternion, ϕc::AbstractVector, 
    timestep)

    X = ∂distance∂x(relative, model.collision, xp, qp, xc, qc)
    Q = ∂distance∂q(relative, model.collision, xp, qp, xc, qc)

    return [X Q]
end

function constraint_jacobian_velocity(relative::Symbol, model::ImpactContact, 
    xp::AbstractVector, vp::AbstractVector, qp::UnitQuaternion, ϕp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::UnitQuaternion, ϕc::AbstractVector, 
    timestep)

    # recover current orientation 
    if relative == :parent
        x = next_position(xp, -vp, timestep)
        ∂x∂v = linear_integrator_jacobian_velocity(x, vp, timestep)
        q = next_orientation(qp, -ϕp, timestep)
        ∂q∂ϕ = rotational_integrator_jacobian_velocity(q, ϕp, timestep)
    elseif relative == :child 
        x = next_position(xc, -vc, timestep)
        ∂x∂v = linear_integrator_jacobian_velocity(x, vc, timestep)
        q = next_orientation(qc, -ϕc, timestep)
        ∂q∂ϕ = rotational_integrator_jacobian_velocity(q, ϕc, timestep)
    end

    # Jacobian
    V = ∂distance∂x(relative, model.collision, xp, qp, xc, qc) * ∂x∂v
    Ω = ∂distance∂q(relative, model.collision, xp, qp, xc, qc) * ∂q∂ϕ

    return [V Ω]
end

function force_mapping(relative::Symbol, model::ImpactContact, 
    xp::AbstractVector, qp::UnitQuaternion, 
    xc::AbstractVector, qc::UnitQuaternion)

    X = contact_normal(relative, model.collision, xp, qp, xc, qc)

    return X
end

function ∂force_mapping∂x(relative::Symbol, jacobian::Symbol,
    model::ImpactContact, 
    xp::AbstractVector, qp::UnitQuaternion, 
    xc::AbstractVector, qc::UnitQuaternion,
    λ::AbstractVector)

    X = ∂contact_normal∂x(relative, jacobian, model.collision, xp, qp, xc, qc, λ)

    return X
end

function ∂force_mapping∂q(relative::Symbol, jacobian::Symbol,
    model::ImpactContact, 
    xp::AbstractVector, qp::UnitQuaternion, 
    xc::AbstractVector, qc::UnitQuaternion,
    λ::AbstractVector)

    X = ∂contact_normal∂q(relative, jacobian, model.collision, xp, qp, xc, qc, λ)

    return X
end



