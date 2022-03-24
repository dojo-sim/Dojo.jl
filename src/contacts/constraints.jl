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

    return impulse_map_jacobian(relative, relative, contact.model,
        get_body(mechanism, contact.parent_id),
        get_body(mechanism, contact.child_id),
        contact.impulses[2],
        mechanism.timestep)
end

function impulses_jacobian_velocity!(mechanism, body::Body, contact::ContactConstraint)
    timestep = mechanism.timestep
    body.state.D -= impulse_map_jacobian_configuration(mechanism, body, contact) * integrator_jacobian_velocity(body, timestep)
    return
end

# off-diagonal terms for linear system
function off_diagonal_jacobians(mechanism, body::Body, contact::RigidContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    Z = szeros(T,N½,6)
    return [Z' -impulse_map(mechanism, contact, body)], [Z; constraint_jacobian_velocity(mechanism, contact, body)]
end

# off-diagonal terms for linear system
function off_diagonal_jacobians(mechanism, contact::RigidContactConstraint{T,N,Nc,Cs,N½}, body::Body) where {T,N,Nc,Cs,N½}
    Z = szeros(T,N½,6)
    return [Z; constraint_jacobian_velocity(mechanism, contact, body)], [Z' -impulse_map(mechanism, contact, body)]
end

# linear system entries
function set_matrix_vector_entries!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry, contact::RigidContactConstraint)
    matrix_entry.value = constraint_jacobian(contact)
    vector_entry.value = vcat(-complementarityμ(mechanism, contact), -constraint(mechanism, contact))
    return
end

# reset variables using cone-specific neutral vector
function reset!(contact::RigidContactConstraint;
    scale=1.0)
    contact.impulses[1] = scale * neutral_vector(contact.model)
    contact.impulses[2] = scale * neutral_vector(contact.model)
    contact.impulses_dual[1] = scale * neutral_vector(contact.model)
    contact.impulses_dual[2] = scale * neutral_vector(contact.model)
    return
end

# dimension of cone
cone_degree(contact::RigidContactConstraint) = sum(cone_degree(contact.model))
