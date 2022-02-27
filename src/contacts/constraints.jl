
# constraints
function constraint_jacobian_configuration(mechanism, contact::ContactConstraint, body::Body)
    return constraint_jacobian_configuration(contact.model, body, nothing, nothing, mechanism.timestep)
end

function constraint_jacobian_velocity(mechanism, contact::ContactConstraint, body::Body)
    return constraint_jacobian_velocity(contact.model, body, nothing, nothing, mechanism.timestep)
end

# impulses
function impulses!(mechanism, body::Body, contact::ContactConstraint)
    body.state.d -= impulse_map(mechanism, contact, body) * contact.impulses[2]
    return
end

function impulse_map(mechanism, contact::ContactConstraint, body::Body)
    return impulse_map(contact.model, body, nothing, nothing, mechanism.timestep)
end

function impulse_map_jacobian_configuration(mechanism, body::Body, contact::ContactConstraint{T}) where T
    x, q = next_configuration(body.state, mechanism.timestep)
    model = contact.model
    X = force_mapping(model, x, q)
    λ = X' * contact.impulses[2]

    ∇Q = skew(model.contact_point - vector_rotate(model.offset, inv(q))) * VRmat(q) * ∂LᵀVᵀmat∂q(λ)
    ∇Q += skew(model.contact_point - vector_rotate(model.offset, inv(q))) * ∂VRmat∂q(LᵀVᵀmat(q) * λ)
    ∇Q += -∂skew∂p(VRmat(q) * LᵀVᵀmat(q) * λ) * ∂rotation_matrix_inv∂q(q, model.offset)
    Z3 = szeros(T,3,3)
    Z4 = szeros(T,3,4)
    return [Z3 Z4;
            Z3 ∇Q]
end

function impulses_jacobian_velocity!(mechanism, body::Body, contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    timestep= mechanism.timestep
    body.state.D -= impulse_map_jacobian_configuration(mechanism, body, contact) * integrator_jacobian_velocity(body, timestep)
    return
end

# off-diagonal terms for linear system
function off_diagonal_jacobians(mechanism, body::Body, contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    Z = szeros(T,N½,6)
    return [Z' -impulse_map(mechanism, contact, body)], [Z; constraint_jacobian_velocity(mechanism, contact, body)]
end

# reset variables using cone-specific neutral vector
function reset!(contact::ContactConstraint{T,N,Nc,Cs,N½}; scale::T=1.0) where {T,N,Nc,Cs,N½}
    contact.impulses_dual[1] = scale * neutral_vector(contact.model)
    contact.impulses_dual[2] = scale * neutral_vector(contact.model)
    contact.impulses[1] = scale * neutral_vector(contact.model)
    contact.impulses[2] = scale * neutral_vector(contact.model)
    return
end

# dimension of cone
cone_degree(contact::ContactConstraint) = sum(cone_degree(contact.model))
