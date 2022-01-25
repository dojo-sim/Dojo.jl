abstract type Contact{T,N} end

getT(bound::Contact{T}) where T = T
Base.length(bound::Contact{T,N}) where {T,N} = N

constraint(bound::Contact, body::Body, λ, Δt) = constraint(bound, next_configuration(body.state, Δt)..., λ)

constraint_jacobian_velocity(bound::Contact, body::Body, id, λ, Δt) = constraint_jacobian_velocity(bound, next_configuration(body.state, Δt)..., current_configuration_velocity(body.state)..., λ, Δt)
constraint_jacobian_configuration(bound::Contact, body::Body, id, λ, Δt) = constraint_jacobian_configuration(bound, next_configuration(body.state, Δt)..., current_configuration_velocity(body.state)..., λ, Δt)

impulse_map(bound::Contact, body::Body, id, λ, Δt) = impulse_map(bound, next_configuration(body.state, Δt)..., λ)

complementarity(mechanism, ineqc::ContactConstraint; scaling=false) = ineqc.γsol[2] .* ineqc.ssol[2]
complementarityμ(mechanism, ineqc::ContactConstraint; scaling=false) = complementarity(mechanism, ineqc, scaling=scaling) - mechanism.μ * neutral_vector(ineqc.constraints[1])

neutral_vector(bound::Contact{T,N}) where {T,N} = sones(T, Int(N/2))

cone_degree(bound::Contact{T,N}) where {T,N} = Int(N/2)
