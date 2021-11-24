abstract type Bound{T,N} end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::Bound{T,N}) where {T,N}
    summary(io, constraint)
end

### General functions
getT(bound::Bound{T}) where T = T
Base.length(bound::Bound{T,N}) where {T,N} = N


### Constraints and derivatives
## Position level constraint wrappers
g(bound::Bound, body::Body, Δt) = g(bound, body.state, Δt)


## Discrete-time position wrappers (for dynamics)
@inline g(bound::Bound, state::State, Δt) = g(bound, posargs3(state, Δt)...)

## Discrete-time position derivatives (for dynamics)
# Wrappers 1
@inline function ∂g∂ʳposa(bound::Bound, body::Body, id, Δt)
    return ∂g∂ʳpos(bound, body.state, Δt)
end

# Wrappers 2
∂g∂ʳpos(bound::Bound, state::State, Δt) = ∂g∂ʳpos(bound, posargs3(state, Δt)...)

# Derivatives accounting for quaternion specialness
@inline function ∂g∂ʳpos(bound::Bound, x::AbstractVector, q::UnitQuaternion)
    X, Q = ∂g∂pos(bound, x, q)
    Q = Q * LVᵀmat(q)
    return [X Q]
end

## Discrete-time velocity derivatives (for dynamics)
# Wrappers 1
@inline function ∂g∂ʳvela(bound::Bound, body::Body, id, Δt)
    return ∂g∂ʳvel(bound, body.state, Δt)
end

# Wrappers 2
∂g∂ʳvel(bound::Bound, state::State, Δt) = ∂g∂ʳvel(bound, posargs3(state, Δt)..., fullargssol(state)..., Δt)

# Derivatives accounting for quaternion specialness
@inline function ∂g∂ʳvel(bound::Bound, x2::AbstractVector, q2::UnitQuaternion,
        x1::AbstractVector, v1::AbstractVector, q1::UnitQuaternion, ω1::AbstractVector, Δt
    )

    X, Q = ∂g∂pos(bound, x2, q2)
    V = X * Δt
    Ω = Q * Lmat(q1) * derivωbar(ω1, Δt) * Δt / 2

    return [V Ω]
end

# signed distance function
function sdf(ineqc::InequalityConstraint{T,N,Nc,Cs}, x::AbstractVector{T}, q::UnitQuaternion{T}) where {T,N,Nc,Cs<:Tuple{<:Bound{T,N}}}
    cont = ineqc.constraints[1]
    return cont.ainv3 * (x + vrotate(cont.p,q) - cont.offset)
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
