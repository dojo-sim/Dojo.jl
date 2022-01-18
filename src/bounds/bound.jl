abstract type Bound{T,N} end

getT(bound::Bound{T}) where T = T
Base.length(bound::Bound{T,N}) where {T,N} = N

g(bound::Bound, body::Body, λ, Δt) = g(bound, posargs3(body.state, Δt)..., λ)

G(bound::Bound, body::Body, id, λ, Δt) = G(bound, posargs3(body.state, Δt)..., λ)

complementarity(mechanism, ineqc::InequalityConstraint; scaling=false) = ineqc.γsol[2] .* ineqc.ssol[2]
complementarityμ(mechanism, ineqc::InequalityConstraint; scaling=false) = complementarity(mechanism, ineqc, scaling=scaling) - mechanism.μ * neutral_vector(ineqc.constraints[1])

neutral_vector(bound::Bound{T,N}) where {T,N} = sones(T, Int(N/2))

cone_degree(bound::Bound{T,N}) where {T,N} = Int(N/2)



