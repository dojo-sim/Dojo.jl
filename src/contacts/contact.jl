abstract type Contact{T,N} end

getT(bound::Contact{T}) where T = T
Base.length(bound::Contact{T,N}) where {T,N} = N

constraint(bound::Contact, body::Body, λ, timestep) = constraint(bound, next_configuration(body.state, timestep)..., λ)

constraint_jacobian_velocity(bound::Contact, body::Body, id, λ, timestep) = constraint_jacobian_velocity(bound, next_configuration(body.state, timestep)..., current_configuration_velocity(body.state)..., λ, timestep)
constraint_jacobian_configuration(bound::Contact, body::Body, id, λ, timestep) = constraint_jacobian_configuration(bound, next_configuration(body.state, timestep)..., current_configuration_velocity(body.state)..., λ, timestep)

impulse_map(bound::Contact, body::Body, id, λ, timestep) = impulse_map(bound, next_configuration(body.state, timestep)..., λ)

complementarity(mechanism, contact::ContactConstraint; scaling=false) = contact.dual[2] .* contact.primal[2]
complementarityμ(mechanism, contact::ContactConstraint; scaling=false) = complementarity(mechanism, contact, scaling=scaling) - mechanism.μ * neutral_vector(contact.constraints[1])

neutral_vector(bound::Contact{T,N}) where {T,N} = sones(T, Int(N/2))

cone_degree(bound::Contact{T,N}) where {T,N} = Int(N/2)
