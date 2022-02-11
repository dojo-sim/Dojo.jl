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
    return constraint_jacobian_configuration(contact.constraints[1], body, nothing, nothing, mechanism.timestep)
end

function impulse_map(mechanism, contact::ContactConstraint, body::Body)
    return impulse_map(contact.model, body, nothing, nothing, mechanism.timestep)
end

@inline function impulses!(mechanism, body::Body, contact::ContactConstraint)
    body.state.d -= impulse_map(mechanism, contact, body) * contact.impulses_dual[2]
    return
end

function impulse_map_jacobian_configuration(mechanism, body::Body, contact::ContactConstraint{T}) where T
    x, q = next_configuration(body.state, mechanism.timestep)
    model = contact.model
    X = force_mapping(model, x, q)
    λ = X' * contact.impulses_dual[2]
    p = model.contact_point

    # Q = skew(p - vrotate(model.offset, inv(q))) * VRmat(q) * LᵀVᵀmat(q) * λ
    ∇Q = skew(p - vrotate(model.offset, inv(q))) * VRmat(q) * ∂qLᵀVᵀmat(λ)
    ∇Q += skew(p - vrotate(model.offset, inv(q))) * ∂qVRmat(LᵀVᵀmat(q) * λ)
    ∇Q += -∂pskew(VRmat(q) * LᵀVᵀmat(q) * λ) * ∂qrotation_matrix_inv(q, model.offset)
    # ∇Q *= LVᵀmat(q)
    Z3 = szeros(T,3,3)
    Z4 = szeros(T,3,4)
    return [Z3 Z4;
            Z3 ∇Q]
end

# @inline function impulse_map_jacobian_velocity(mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½}, body::Body) where {T,N,Nc,Cs,N½}
#     timestep = mechanism.timestep
#     x3, q3 = next_configuration(body.state, timestep)
#     x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)
#     D = szeros(T,6,6)
#     for i=1:Nc
#         bnd = contact.constraints[i]
#         p = bnd.p
#         offset = bnd.offset
#
#         X = force_mapping(bnd)
#         λ = X' * contact.dual[2]
#
#         ∇ = ∂pskew(VRmat(q3) * LᵀVᵀmat(q3) * λ) * -∂vrotate∂q(offset, inv(q3)) * Tmat()
#         ∇ += skew(p - vrotate(offset, inv(q3))) * ∂qVRmat(LᵀVᵀmat(q3) * λ)
#         ∇ += skew(p - vrotate(offset, inv(q3))) * VRmat(q3) * ∂qLᵀVᵀmat(λ)
#         ∇ *= rotational_integrator_jacobian_velocity(q2, ϕ25, timestep)
#         D = [szeros(T,6,3) [szeros(T,3,3); ∇]]
#     end
#     return D
# end
@inline function impulses_jacobian_velocity!(mechanism, body::Body, contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    timestep = mechanism.timestep
    body.state.D -= impulse_map_jacobian_configuration(mechanism, body, contact) *
        integrator_jacobian_velocity(body, timestep)
    # x3, q3 = next_configuration(body.state, timestep)
    # x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)
    #
    # for i=1:Nc
    #     model = contact.model
    #     contact_point = model.contact_point
    #     offset = model.offset
    #
    #     X = force_mapping(model)
    #     λ = X' * contact.impulses_dual[2]
    #
    #     ∇ = ∂pskew(VRmat(q3) * LᵀVᵀmat(q3) * λ) * -∂vrotate∂q(offset, inv(q3)) * Tmat()
    #     ∇ += skew(contact_point - vrotate(offset, inv(q3))) * ∂qVRmat(LᵀVᵀmat(q3) * λ)
    #     ∇ += skew(contact_point - vrotate(offset, inv(q3))) * VRmat(q3) * ∂qLᵀVᵀmat(λ)
    #     ∇ *= rotational_integrator_jacobian_velocity(q2, ϕ25, timestep)
    #     body.state.D -= [szeros(T,6,3) [szeros(T,3,3); ∇]]
    # end
    return
end

@inline function off_diagonal_jacobians(mechanism, body::Body, contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    Z = szeros(T,N½,6)
    return [Z' -impulse_map(mechanism, contact, body)], [Z; constraint_jacobian_velocity(mechanism, contact, body)]
end

function reset!(contact::ContactConstraint{T,N,Nc,Cs,N½}; scale::T=1.0) where {T,N,Nc,Cs,N½}
    contact.impulses[1] = scale * neutral_vector(contact.model)
    contact.impulses[2] = scale * neutral_vector(contact.model)
    contact.impulses_dual[1] = scale * neutral_vector(contact.model)
    contact.impulses_dual[2] = scale * neutral_vector(contact.model)
    return
end

cone_degree(contact::ContactConstraint) = sum(cone_degree(contact.model))
