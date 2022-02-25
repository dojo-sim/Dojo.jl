abstract type Contact{T,N} end

# constraints
constraint(model::Contact, body::Body, λ, timestep) = constraint(model, next_configuration(body.state, timestep)..., λ)

# constraint Jacobians
constraint_jacobian_configuration(model::Contact, body::Body, id, λ, timestep) = constraint_jacobian_configuration(model, next_configuration(body.state, timestep)..., current_configuration_velocity(body.state)..., λ, timestep)
constraint_jacobian_velocity(model::Contact, body::Body, id, λ, timestep) = constraint_jacobian_velocity(model, next_configuration(body.state, timestep)..., current_configuration_velocity(body.state)..., λ, timestep)

# impulses
impulse_map(model::Contact, body::Body, id, λ, timestep) = impulse_map(model, next_configuration(body.state, timestep)..., λ)

function impulse_map(model::Contact, x::AbstractVector, q::UnitQuaternion, λ)
    X = force_mapping(model, x, q)
    Q = - X * q * skew(model.contact_point - vector_rotate(model.offset, inv(q)))
    return [X'; Q']
end

# force mapping 
function force_mapping(model::Contact, x::AbstractVector, q::UnitQuaternion)
    X = [model.surface_normal_projector;
         szeros(1,3);
         model.surface_projector]
    return X
end

# utilities
Base.length(model::Contact{T,N}) where {T,N} = N
neutral_vector(model::Contact{T,N}) where {T,N} = sones(T, Int(N/2))
cone_degree(model::Contact{T,N}) where {T,N} = Int(N/2)




