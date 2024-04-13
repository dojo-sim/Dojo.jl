"""
    ImpactContact{T,N} <: Contact{T,N}

    contact object for impact (i.e., no friction)

    collision: Collision
"""
mutable struct ImpactContact{T,N} <: Contact{T,N}
    friction_parameterization::SMatrix{0,2,T,0}
    collision::Collision{T,0,3,0}
end

# function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, contact::ImpactContact)
#     summary(io, contact)
#     println(io, "")
#     println(io, "friction_parameterization: "*string(contact.friction_parameterization))
#     println(io, "collision:                 "*string(contact.collision))
# end

function ImpactContact(body::Body{T}, normal::AbstractVector; 
    contact_origin=szeros(T, 3), 
    contact_radius=0.0) where T
    
    # contact directions
    V1, V2, V3 = orthogonal_columns(normal) #
    A = [V1 V2 V3]
    Ainv = inv(A)
    contact_normal = Ainv[3, SA[1; 2; 3]]'

    # friction parametrization
    parameterization = szeros(T, 0, 2)

    # collision 
    collision = SphereHalfSpaceCollision(szeros(T, 0, 3), contact_normal, SVector{3}(contact_origin), contact_radius)

    ImpactContact{Float64,2}(parameterization, collision), body.id, 0
end

function ImpactContact(body::Body{T}, normals::AbstractVector{<:AbstractVector}; 
    contact_origins::AbstractVector=[szeros(T, 3) for i=1:length(normals)],
    contact_radii::AbstractVector=[0.0 for i=1:length(normals)]) where T

    @assert length(normals) == length(contact_origins) == length(contact_radii)
    [ImpactContact(body, normals[i]; contact_origin=contact_origins[i], contact_radius=contact_radii[i]) for i in eachindex(normals)]
end

function ImpactContact(bodies::AbstractVector{Body{T}}, normals::AbstractVector{<:AbstractVector}; 
    contact_origins::AbstractVector=[szeros(T, 3) for i=1:length(normals)],
    contact_radii::AbstractVector=[0.0 for i=1:length(normals)]) where T

    @assert length(bodies) == length(normals) == length(contact_origins) == length(contact_radii)
    [ImpactContact(bodies[i], normals[i]; contact_origin=contact_origins[i], contact_radius=contact_radii[i]) for i in eachindex(bodies)]
end

function ImpactContact(body1::Body{T}, body2::Body{T}) where T
    # friction parametrization
    parameterization = szeros(T, 0, 2)

    origin_parent = szeros(3)
    origin_child = szeros(3)
    radius_parent = 1
    radius_child = 1
    primitive1 = body1.shape.primitive
    primitive2 = body2.shape.primitive
    α = 0
    intersection_point = szeros(3)

    collision = GeneralCollision{T,0,3,0}(origin_parent,origin_child,radius_parent,radius_child,primitive1,primitive2,α,intersection_point)

    ImpactContact{Float64,2}(parameterization, collision), body1.id, body2.id
end

function ImpactContact(bodies::AbstractVector{Body{T}}) where T
    impacts = Tuple[]
    for i=1:length(bodies)
        for j=i+1:length(bodies)
            impacts = [impacts;ImpactContact(bodies[i],bodies[j])]
        end
    end
    
    return impacts
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
    xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ωp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ωc::AbstractVector, 
    timestep)

    X = ∂distance∂x(relative, model.collision, xp, qp, xc, qc)
    Q = ∂distance∂q(relative, model.collision, xp, qp, xc, qc)

    return [X Q]
end

function constraint_jacobian_velocity(relative::Symbol, model::ImpactContact, 
    xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ωp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ωc::AbstractVector, 
    timestep)

    # recover current orientation 
    if relative == :parent
        x = next_position(xp, -vp, timestep)
        ∂x∂v = linear_integrator_jacobian_velocity(x, vp, timestep)
        q = next_orientation(qp, -ωp, timestep)
        ∂q∂ω = rotational_integrator_jacobian_velocity(q, ωp, timestep)
    elseif relative == :child 
        x = next_position(xc, -vc, timestep)
        ∂x∂v = linear_integrator_jacobian_velocity(x, vc, timestep)
        q = next_orientation(qc, -ωc, timestep)
        ∂q∂ω = rotational_integrator_jacobian_velocity(q, ωc, timestep)
    end

    # Jacobian
    V = ∂distance∂x(relative, model.collision, xp, qp, xc, qc) * ∂x∂v
    Ω = ∂distance∂q(relative, model.collision, xp, qp, xc, qc) * ∂q∂ω

    return [V Ω]
end

function force_mapping(relative::Symbol, model::ImpactContact, 
    xp::AbstractVector, qp::Quaternion, 
    xc::AbstractVector, qc::Quaternion)

    X = contact_normal(model.collision, xp, qp, xc, qc)'

    if relative == :parent 
        return X 
    elseif relative == :child 
        return -1.0 * X 
    end
end

function ∂force_mapping_jvp∂x(relative::Symbol, jacobian::Symbol,
    model::ImpactContact, 
    xp::AbstractVector, qp::Quaternion, 
    xc::AbstractVector, qc::Quaternion,
    λ::AbstractVector)

    X = ∂contact_normal_transpose∂x(jacobian, model.collision, xp, qp, xc, qc) * λ[1]

    if relative == :parent 
        return X 
    elseif relative == :child 
        return -1.0 * X 
    end
end

function ∂force_mapping_jvp∂q(relative::Symbol, jacobian::Symbol,
    model::ImpactContact, 
    xp::AbstractVector, qp::Quaternion, 
    xc::AbstractVector, qc::Quaternion,
    λ::AbstractVector)

    X = ∂contact_normal_transpose∂q(jacobian, model.collision, xp, qp, xc, qc) * λ[1]

    if relative == :parent 
        return X 
    elseif relative == :child 
        return -1.0 * X 
    end
end



