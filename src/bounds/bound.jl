abstract type Bound{T,N} end

### General functions
getT(bound::Bound{T}) where T = T
Base.length(bound::Bound{T,N}) where {T,N} = N


## Position level constraint wrappers
g(bound::Bound, body::Body, λ, Δt) = g(bound, posargs3(body.state, Δt)..., λ)

## Discrete-time position derivatives (for dynamics)
# Wrappers 1
@inline function Ga(bound::Bound, body::Body, id, λ, Δt)
    return G(bound, posargs3(body.state, Δt)..., λ)
end

# # Derivatives accounting for quaternion specialness
# @inline function G(bound::Bound, x::AbstractVector, q::UnitQuaternion, λ)
#     X, Q = ∂g∂pos(bound, x, q, λ)
#     Q = Q * LVᵀmat(q)
#     return [X Q]
# end

# ## Discrete-time velocity derivatives (for dynamics)
# # Wrappers 1
# @inline function ∂g∂ʳvela(bound::Bound, body::Body, id, λ, Δt)
#     return ∂g∂ʳvel(bound, body.state, λ, Δt)
# end

# # Wrappers 2
# ∂g∂ʳvel(bound::Bound, state::State, λ, Δt) = ∂g∂ʳvel(bound, posargs3(state, Δt)..., fullargssol(state)..., λ, Δt)

# # Derivatives accounting for quaternion specialness
# @inline function ∂g∂ʳvel(bound::Bound, x2::AbstractVector, q2::UnitQuaternion,
#         x1::AbstractVector, v1::AbstractVector, q1::UnitQuaternion, ω1::AbstractVector, λ, Δt
#     )

#     X, Q = ∂g∂pos(bound, x2, q2, λ)
#     V = X * Δt
#     Ω = Q * Lmat(q1) * derivωbar(ω1, Δt) * Δt / 2

#     return [V Ω]
# end

# signed distance function
function sdf(ineqc::InequalityConstraint{T,N,Nc,Cs}, x::AbstractVector{T}, q::UnitQuaternion{T}) where {T,N,Nc,Cs<:Tuple{<:Bound{T,N}}}
    cont = ineqc.constraints[1]
    return cont.ainv3 * (x + vrotate(cont.p, q) - cont.offset)
end

# contact location
function contact_location(mechanism::Mechanism)
    return [contact_location(mech, ineqc) for ineqc in mechanism.ineqconstraints]
end

function contact_location(mechanism::Mechanism, ineqc::InequalityConstraint)
    bodies = collect(mechanism.bodies)
    body = bodies[findfirst(x -> x.id == ineqc.parentid, bodies)]
    return contact_location(ineqc, body)
end

function contact_location(ineqc::InequalityConstraint{T,N,Nc,Cs}, body::Body) where {T,N,Nc,Cs<:Tuple{<:Bound{T,N}}}
    x = body.state.x2[1]
    q = body.state.q2[1]
    return contact_location(ineqc, x, q)
end

function contact_location(ineqc::InequalityConstraint{T,N,Nc,Cs}, x::AbstractVector{T}, q::UnitQuaternion{T}) where {T,N,Nc,Cs<:Tuple{<:Bound{T,N}}}
    cont = ineqc.constraints[1]
    return x + vrotate(cont.p,q) - cont.offset
end

function get_sdf(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, storage::Storage) where {T,Nn,Ne,Nb,Ni}
    N = length(storage.x[1])
    d = []
    for ineqc in mechanism.ineqconstraints
        ibody = getbody(mechanism, ineqc.parentid).id - Ne
        push!(d, [sdf(ineqc, storage.x[ibody][i], storage.q[ibody][i]) for i = 1:length(storage.x[1])])
    end
    return d
end

