# # constraint Jacobian
# function constraint_jacobian_configuration(mechanism, contact::SoftContactConstraint, body::Body)
#     relative = (body.id == contact.parent_id ? :parent : :child)
#
#     timestep = mechanism.timestep
#     pbody = get_body(mechanism, contact.parent_id)
#     cbody = get_body(mechanism, contact.child_id)
#
#     return constraint_jacobian_configuration(relative,
#         contact.model,
#         next_configuration_velocity(pbody.state, timestep)...,
#         next_configuration_velocity(cbody.state, timestep)...,
#         mechanism.timestep)
# end
#
# function constraint_jacobian_velocity(mechanism, contact::SoftContactConstraint, body::Body)
#     relative = (body.id == contact.parent_id ? :parent : :child)
#
#     timestep = mechanism.timestep
#     pbody = get_body(mechanism, contact.parent_id)
#     cbody = get_body(mechanism, contact.child_id)
#
#     return constraint_jacobian_velocity(relative,
#         contact.model,
#         next_configuration_velocity(pbody.state, timestep)...,
#         next_configuration_velocity(cbody.state, timestep)...,
#         mechanism.timestep)
# end
#
# # impulses
# function impulses!(mechanism, body::Body, contact::SoftContactConstraint)
#     body.state.d -= impulse_map(mechanism, contact, body) * contact.impulses[2]
#     return
# end
#
# function impulse_map(mechanism, contact::SoftContactConstraint, body::Body)
#     relative = (body.id == contact.parent_id ? :parent : :child)
#     pbody = get_body(mechanism, contact.parent_id)
#     cbody = get_body(mechanism, contact.child_id)
#     return impulse_map(relative, contact.model, pbody, cbody, mechanism.timestep)
# end
#
# function impulse_map_jacobian_configuration(mechanism, body::Body{T}, contact::SoftContactConstraint{T}) where T
#     relative = (body.id == contact.parent_id ? :parent : :child)
#
#     return impulse_map_jacobian(relative, relative, contact.model,
#         get_body(mechanism, contact.parent_id),
#         get_body(mechanism, contact.child_id),
#         contact.impulses[2],
#         mechanism.timestep)
# end
#
# function impulses_jacobian_velocity!(mechanism, body::Body, contact::SoftContactConstraint)
#     timestep = mechanism.timestep
#     @warn "not sure about this one"
#     body.state.D -= impulse_map_jacobian_configuration(mechanism, body, contact) * integrator_jacobian_velocity(body, timestep)
#     return
# end

# off-diagonal terms for linear system
function off_diagonal_jacobians(mechanism, body::Body, contact::SoftContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs}
    return -impulse_map(mechanism, contact, body), constraint_jacobian_velocity(mechanism, contact, body)
end

# off-diagonal terms for linear system
function off_diagonal_jacobians(mechanism, contact::SoftContactConstraint{T,N,Nc,Cs}, body::Body) where {T,N,Nc,Cs}
    return constraint_jacobian_velocity(mechanism, contact, body), -impulse_map(mechanism, contact, body)
end

# linear system entries
function set_matrix_vector_entries!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry, contact::SoftContactConstraint)
    matrix_entry.value = constraint_jacobian(contact)
    vector_entry.value = -constraint(mechanism, contact)
    return
end

# reset variables using cone-specific neutral vector
function reset!(contact::SoftContactConstraint{T,N}; scale=0.0) where {T,N}
    contact.impulses[1] = szeros(T,N)
    contact.impulses[2] = szeros(T,N)
    return
end

# initialization
function initialize!(contact::SoftContactConstraint)
    return nothing
end


# complementarity
complementarity(mechanism, contact::SoftContactConstraint{T,N}; scaling=false) where {T,N} = szeros(T,N)
complementarityμ(mechanism, contact::SoftContactConstraint{T,N}; scaling=false) where {T,N} = szeros(T,N)


# cone line search
function cone_line_search!(α, mechanism, contact::SoftContactConstraint,
        vector_entry::Entry, τort, τsoc; scaling::Bool=false)
    return α
end


# centering
function centering!(ν, νaff, n, mechanism, contact::SoftContactConstraint, vector_entry::Entry, αaff)
    return ν, νaff, n
end


# candidate step
function candidate_step!(α::T, contact::SoftContactConstraint{T}, vector_entry::Entry, scale) where T
    contact.impulses[2] = contact.impulses[1] + 1 / (2^scale) * α * vector_entry.value
    return
end


# update
function update!(contact::SoftContactConstraint)
    contact.impulses[1] = contact.impulses[2]
    return
end
