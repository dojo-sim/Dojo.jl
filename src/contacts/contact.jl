abstract type Contact{T,N} end

getT(bound::Contact{T}) where T = T
Base.length(bound::Contact{T,N}) where {T,N} = N

g(bound::Contact, body::Body, λ, Δt) = g(bound, next_configuration(body.state, Δt)..., λ)

∂g∂v(bound::Contact, body::Body, id, λ, Δt) = ∂g∂v(bound, next_configuration(body.state, Δt)..., current_configuration_velocity(body.state)..., λ, Δt)
∂g∂z(bound::Contact, body::Body, id, λ, Δt) = ∂g∂z(bound, next_configuration(body.state, Δt)..., current_configuration_velocity(body.state)..., λ, Δt)

G(bound::Contact, body::Body, id, λ, Δt) = G(bound, next_configuration(body.state, Δt)..., λ)

complementarity(mechanism, ineqc::ContactConstraint; scaling=false) = ineqc.γsol[2] .* ineqc.ssol[2]
complementarityμ(mechanism, ineqc::ContactConstraint; scaling=false) = complementarity(mechanism, ineqc, scaling=scaling) - mechanism.μ * neutral_vector(ineqc.constraints[1])

neutral_vector(bound::Contact{T,N}) where {T,N} = sones(T, Int(N/2))

cone_degree(bound::Contact{T,N}) where {T,N} = Int(N/2)
