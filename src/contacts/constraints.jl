# constraint Jacobian 
function constraint_jacobian_configuration(mechanism, contact::ContactConstraint, body::Body)
    relative = (body.id == contact.parent_id ? :parent : :child)
    
    timestep = mechanism.timestep
    pbody = get_body(mechanism, contact.parent_id) 
    cbody = get_body(mechanism, contact.child_id)
    
    return constraint_jacobian_configuration(relative,
        contact.model, 
        next_configuration_velocity(pbody.state, timestep)..., 
        next_configuration_velocity(cbody.state, timestep)..., 
        mechanism.timestep)
end

function constraint_jacobian_velocity(mechanism, contact::ContactConstraint, body::Body)
    relative = (body.id == contact.parent_id ? :parent : :child)

    timestep = mechanism.timestep
    pbody = get_body(mechanism, contact.parent_id) 
    cbody = get_body(mechanism, contact.child_id)

    return constraint_jacobian_velocity(relative,
        contact.model, 
        next_configuration_velocity(pbody.state, timestep)..., 
        next_configuration_velocity(cbody.state, timestep)..., 
        mechanism.timestep)
end

# impulses
function impulses!(mechanism, body::Body, contact::ContactConstraint)
    body.state.d -= impulse_map(mechanism, contact, body) * contact.impulses[2]
    return
end

function impulse_map(mechanism, contact::ContactConstraint, body::Body)
    relative = (body.id == contact.parent_id ? :parent : :child)
    pbody = get_body(mechanism, contact.parent_id) 
    cbody = get_body(mechanism, contact.child_id)
    return impulse_map(relative, contact.model, pbody, cbody, mechanism.timestep)
end

function impulse_map_jacobian_configuration(mechanism, body::Body{T}, contact::ContactConstraint{T}) where T
    relative = (body.id == contact.parent_id ? :parent : :child)

    # contact model 
    model = contact.model

    # parent 
    xp, qp = next_configuration(get_body(mechanism, contact.parent_id).state, mechanism.timestep)

    # child
    xc, qc = next_configuration(get_body(mechanism, contact.child_id).state, mechanism.timestep)

    # contact impulse
    X = force_mapping(relative, model, xp, qp, xc, qc)
    λ = X' * contact.impulses[2]

    # offset 
    offset = model.collision.contact_normal' * model.collision.contact_radius

    # Jacobian 
    Z3 = szeros(T,3,3)
    Z4 = szeros(T,3,4)

    ∇Q = skew(model.collision.contact_origin - vector_rotate(offset, inv(qp))) * VRmat(qp) * ∂LᵀVᵀmat∂q(λ)
    ∇Q += skew(model.collision.contact_origin - vector_rotate(offset, inv(qp))) * ∂VRmat∂q(LᵀVᵀmat(qp) * λ)
    ∇Q += -∂skew∂p(VRmat(qp) * LᵀVᵀmat(qp) * λ) * ∂rotation_matrix_inv∂q(qp, offset)
    
    return [Z3 Z4;
            Z3 ∇Q]
end

function impulses_jacobian_velocity!(mechanism, body::Body, contact::ContactConstraint)
    timestep = mechanism.timestep
    body.state.D -= impulse_map_jacobian_configuration(mechanism, body, contact) * integrator_jacobian_velocity(body, timestep)
    return
end

# off-diagonal terms for linear system
function off_diagonal_jacobians(mechanism, body::Body, contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    Z = szeros(T,N½,6)
    return [Z' -impulse_map(mechanism, contact, body)], [Z; constraint_jacobian_velocity(mechanism, contact, body)]
end

# off-diagonal terms for linear system
function off_diagonal_jacobians(mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½}, body::Body) where {T,N,Nc,Cs,N½}
    Z = szeros(T,N½,6)
    return [Z; constraint_jacobian_velocity(mechanism, contact, body)], [Z' -impulse_map(mechanism, contact, body)]
end

# linear system entries
function set_matrix_vector_entries!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry, contact::ContactConstraint)
    matrix_entry.value = constraint_jacobian(contact)
    vector_entry.value = vcat(-complementarityμ(mechanism, contact), -constraint(mechanism, contact))
    return
end

# reset variables using cone-specific neutral vector
function reset!(contact::ContactConstraint; 
    scale=1.0)
    contact.impulses[1] = scale * neutral_vector(contact.model)
    contact.impulses[2] = scale * neutral_vector(contact.model)
    contact.impulses_dual[1] = scale * neutral_vector(contact.model)
    contact.impulses_dual[2] = scale * neutral_vector(contact.model)
    return
end

# dimension of cone
cone_degree(contact::ContactConstraint) = sum(cone_degree(contact.model))
