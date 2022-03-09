"""
    ImpactContact{T,N} <: Contact{T,N}

    contact object for impact (i.e., no friction)

    surface_normal_projector: mapping from world frame to surface tangent frame 
    contact_point: position of contact on Body relative to center of mass 
    offset: position of contact relative to contact_point
"""
mutable struct ImpactContact{T,N} <: Contact{T,N}
    surface_normal_projector::Adjoint{T,SVector{3,T}}
    contact_point::SVector{3,T}
    offset::SVector{3,T}

    function ImpactContact(body::Body{T}, normal::AbstractVector; 
        contact_point=szeros(T, 3), 
        offset::AbstractVector=szeros(T, 3)) where T
        V1, V2, V3 = orthogonal_columns(normal) #
        A = [V1 V2 V3]
        Ainv = inv(A)
        surface_normal_projector = Ainv[3,SA[1; 2; 3]]'
        new{Float64,2}(surface_normal_projector, contact_point, offset)
    end
end

function constraint(mechanism, contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:ImpactContact{T,N}}
    # contact model
    model = contact.model

    # parent
    pbody = get_body(mechanism, contact.parent_id)
    xp3, qp3 = next_configuration(pbody.state, mechanism.timestep)

    # child
    cbody = get_body(mechanism, contact.child_id)
    xc3, qc3 = next_configuration(cbody.state, mechanism.timestep)

    # constraint 
    SVector{1,T}(distance(model, xp3, qp3, xc3, qc3) - contact.impulses_dual[2][1])
end

function constraint_jacobian_configuration(model::ImpactContact, 
    xp3::AbstractVector, vp25::AbstractVector, qp3::UnitQuaternion, ϕp25::AbstractVector,
    xc3::AbstractVector, vc25::AbstractVector, qc3::UnitQuaternion, ϕc25::AbstractVector, 
    timestep)

    X = ∂distance∂xp(model, xp3, qp3, x3c, qc3)
    Q = ∂distance∂qp(model, xp3, qp3, x3c, qc3)

    return [X Q]
end

function constraint_jacobian_velocity(model::ImpactContact, 
    xp3::AbstractVector, vp25::AbstractVector, qp3::UnitQuaternion, ϕp25::AbstractVector,
    xc3::AbstractVector, vc25::AbstractVector, qc3::UnitQuaternion, ϕc25::AbstractVector, 
    timestep)

    # recover current orientation 
    qp2 = next_orientation(qp3, -ϕp25, timestep)

    # Jacobian
    V = ∂distance∂xp(model, xp3, qp3, xc3, qc3) * timestep
    Ω = ∂distance∂qp(model, xp3, qp3, xc3, qc3) * rotational_integrator_jacobian_velocity(qp2, ϕp25, timestep)

    return [V Ω]
end

function force_mapping(model::ImpactContact, xp::AbstractVector, qp::UnitQuaternion, xc::AbstractVector, qc::UnitQuaternion)
    X = model.surface_normal_projector
    return X
end

function complementarity_jacobian(contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:ImpactContact{T,N},N½}
    γ = contact.impulses[2] + REG * neutral_vector(contact.model)
    s = contact.impulses_dual[2] + REG * neutral_vector(contact.model)
    ∇s = hcat(γ, -Diagonal(sones(N½)))
    ∇γ = hcat(s, -Diagonal(szeros(N½)))
    return [∇s ∇γ]
end

