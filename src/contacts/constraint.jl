mutable struct ContactConstraint{T,N,Nc,Cs,N½} <: Constraint{T,N}
    # ID
    id::Int64
    name::Symbol

    # contact model
    model::Cs

    # neighbor IDs
    parent_id::Int
    child_id::Int

    # variables
    impulses::Vector{SVector{N½,T}}
    impulses_dual::Vector{SVector{N½,T}}

    function ContactConstraint(data; name::Symbol=Symbol("contact_" * randstring(4)))
        model, parent_id, _ = data
        T = getT(model)

        N = length(model)
        N½ = Int64(N/2)

        impulses = [neutral_vector(model) for i = 1:2]
        impulses_dual = [neutral_vector(model) for i = 1:2]
        new{T,N,1,typeof(model),N½}(getGlobalID(), name, model, parent_id, 0, impulses, impulses_dual)
    end
end

function constraint_jacobian_velocity(mechanism, contact::ContactConstraint, body::Body)
    return constraint_jacobian_velocity(contact.model, body, nothing, nothing, mechanism.timestep)
end

function constraint_jacobian_configuration(mechanism, contact::ContactConstraint, body::Body)
    return constraint_jacobian_configuration(contact.model, body, nothing, nothing, mechanism.timestep)
end

function impulse_map(mechanism, contact::ContactConstraint, body::Body)
    return impulse_map(contact.model, body, nothing, nothing, mechanism.timestep)
end

@inline function impulses!(mechanism, body::Body, contact::ContactConstraint)
    body.state.d -= impulse_map(mechanism, contact, body) * contact.impulses[2]
    return
end

function impulse_map_jacobian_configuration(mechanism, body::Body, contact::ContactConstraint{T}) where T
    x, q = next_configuration(body.state, mechanism.timestep)
    model = contact.model
    X = force_mapping(model, x, q)
    λ = X' * contact.impulses[2]

    # Q = skew(model.contact_point - vrotate(model.offset, inv(q))) * VRmat(q) * LᵀVᵀmat(q) * λ
    ∇Q = skew(model.contact_point - vrotate(model.offset, inv(q))) * VRmat(q) * ∂qLᵀVᵀmat(λ)
    ∇Q += skew(model.contact_point - vrotate(model.offset, inv(q))) * ∂qVRmat(LᵀVᵀmat(q) * λ)
    ∇Q += -∂pskew(VRmat(q) * LᵀVᵀmat(q) * λ) * ∂qrotation_matrix_inv(q, model.offset)
    # ∇Q *= LVᵀmat(q)
    Z3 = szeros(T,3,3)
    Z4 = szeros(T,3,4)
    return [Z3 Z4;
            Z3 ∇Q]
end

@inline function impulses_jacobian_velocity!(mechanism, body::Body, contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    timestep = mechanism.timestep
    body.state.D -= impulse_map_jacobian_configuration(mechanism, body, contact) *
        integrator_jacobian_velocity(body, timestep)
    return
end

@inline function off_diagonal_jacobians(mechanism, body::Body, contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    Z = szeros(T,N½,6)
    return [Z' -impulse_map(mechanism, contact, body)], [Z; constraint_jacobian_velocity(mechanism, contact, body)]
end

function reset!(contact::ContactConstraint{T,N,Nc,Cs,N½}; scale::T=1.0) where {T,N,Nc,Cs,N½}
    contact.impulses_dual[1] = scale * neutral_vector(contact.model)
    contact.impulses_dual[2] = scale * neutral_vector(contact.model)
    contact.impulses[1] = scale * neutral_vector(contact.model)
    contact.impulses[2] = scale * neutral_vector(contact.model)
    return
end

cone_degree(contact::ContactConstraint) = sum(cone_degree(contact.model))
