mutable struct ContactConstraint{T,N,Nc,Cs,N½} <: Constraint{T,N}
    # ID
    id::Int64
    name::Symbol

    # currently only single constraint and child
    constraints::Cs

    # neighbor IDs
    parent_id::Int
    child_id::Int

    # variables
    primal::Vector{SVector{N½,T}} # holds the slack variable
    dual::Vector{SVector{N½,T}} # holds the dual of the slack variable

    function ContactConstraint(data; name::Symbol=Symbol("contact_" * randstring(4)))
        bound, parent_id, _ = data
        T = getT(bound)

        constraint = Tuple([bound])
        N = length(constraint[1])
        N½ = Int64(N/2)

        primal = [neutral_vector(bound) for i = 1:2]
        dual = [neutral_vector(bound) for i = 1:2]
        new{T,N,1,typeof(constraint),N½}(getGlobalID(), name, constraint, parent_id, 0, primal, dual)
    end
end

function constraint_jacobian_velocity(mechanism, contact::ContactConstraint, body::Body)
    return constraint_jacobian_velocity(contact.constraints[1], body, nothing, nothing, mechanism.timestep)
end

function constraint_jacobian_configuration(mechanism, contact::ContactConstraint, body::Body)
    return constraint_jacobian_configuration(contact.constraints[1], body, nothing, nothing, mechanism.timestep)
end

function impulse_map(mechanism, contact::ContactConstraint, body::Body)
    return impulse_map(contact.constraints[1], body, nothing, nothing, mechanism.timestep)
end

@inline function impulses!(mechanism, body::Body, contact::ContactConstraint)
    body.state.d -= impulse_map(mechanism, contact, body)' * contact.dual[2]
    return
end

@inline function impulses_jacobian_velocity!(mechanism, body::Body, contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    timestep = mechanism.timestep
    x3, q3 = next_configuration(body.state, timestep)
    x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)

    for i=1:Nc
        bnd = contact.constraints[i]
        p = bnd.p
        offset = bnd.offset

        X = force_mapping(bnd)
        λ = X' * contact.dual[2]

        ∇ = ∂pskew(VRmat(q3) * LᵀVᵀmat(q3) * λ) * -∂vrotate∂q(offset, inv(q3)) * Tmat()
        ∇ += skew(p - vrotate(offset, inv(q3))) * ∂qVRmat(LᵀVᵀmat(q3) * λ)
        ∇ += skew(p - vrotate(offset, inv(q3))) * VRmat(q3) * ∂qLᵀVᵀmat(λ)
        ∇ *= rotational_integrator_jacobian_velocity(q2, ϕ25, timestep)
        body.state.D -= [szeros(T,6,3) [szeros(T,3,3); ∇]]
    end
    return
end

@inline function off_diagonal_jacobians(mechanism, body::Body, contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    Z = szeros(T,N½,6)
    return [Z; -impulse_map(mechanism, contact, body)]', [Z; constraint_jacobian_velocity(mechanism, contact, body)]
end

function reset!(contact::ContactConstraint{T,N,Nc,Cs,N½}; scale::T=1.0) where {T,N,Nc,Cs,N½}
    contact.primal[1] = scale * neutral_vector(contact.constraints[1])
    contact.primal[2] = scale * neutral_vector(contact.constraints[1])
    contact.dual[1] = scale * neutral_vector(contact.constraints[1])
    contact.dual[2] = scale * neutral_vector(contact.constraints[1])
    return
end

cone_degree(contact::ContactConstraint) = sum(cone_degree.(contact.constraints))
