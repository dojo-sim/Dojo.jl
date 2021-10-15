abstract type AbstractJoint{T,N} end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::AbstractJoint{T,N}) where {T,N}
    summary(io, constraint)
end

### General functions
getT(joint::AbstractJoint{T}) where T = T
Base.length(joint::AbstractJoint{T,N}) where {T,N} = N
Base.zero(joint::AbstractJoint{T,N}) where {T,N} = szeros(T, N, 6)


### Constraints and derivatives
## Position level constraint wrappers
@inline g(joint::AbstractJoint, body1::Body, body2::Body, Δt) = constraintmat(joint) * g(joint, body1.state, body2.state, Δt)
@inline g(joint::AbstractJoint, body1::Origin, body2::Body, Δt) = constraintmat(joint) * g(joint, body2.state, Δt)
@inline g(joint::AbstractJoint, body1::Body, body2::Body) = constraintmat(joint) * g(joint, body1.state, body2.state)
@inline g(joint::AbstractJoint, body1::Origin, body2::Body) = constraintmat(joint) * g(joint, body2.state)

## Discrete-time position derivatives (for dynamics)
# Wrappers 1
@inline function ∂g∂ʳposa(joint::AbstractJoint, body1::Body, body2::Body, childid)
    if body2.id == childid
        return constraintmat(joint) * ∂g∂ʳposa(joint, body1.state, body2.state)
    else
        return zero(joint)
    end
end
@inline function ∂g∂ʳposb(joint::AbstractJoint, body1::Body, body2::Body, childid)
    if body2.id == childid
        return constraintmat(joint) * ∂g∂ʳposb(joint, body1.state, body2.state)
    else
        return zero(joint)
    end
end
@inline function ∂g∂ʳposb(joint::AbstractJoint, body1::Origin, body2::Body, childid)
    if body2.id == childid
        return constraintmat(joint) * ∂g∂ʳposb(joint, body2.state)
    else
        return zero(joint)
    end
end

# Wrappers 2
∂g∂ʳposa(joint::AbstractJoint, statea::State, stateb::State) = ∂g∂ʳposa(joint, posargsk(statea)..., posargsk(stateb)...)
∂g∂ʳposb(joint::AbstractJoint, statea::State, stateb::State) = ∂g∂ʳposb(joint, posargsk(statea)..., posargsk(stateb)...)
∂g∂ʳposb(joint::AbstractJoint, stateb::State) = ∂g∂ʳposb(joint, posargsk(stateb)...)

# Derivatives accounting for quaternion specialness
@inline function ∂g∂ʳposa(joint::AbstractJoint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    X, Q = ∂g∂posa(joint, xa, qa, xb, qb)
    Q = Q * LVᵀmat(qa)

    return [X Q]
end
@inline function ∂g∂ʳposb(joint::AbstractJoint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    X, Q = ∂g∂posb(joint, xa, qa, xb, qb)
    Q = Q * LVᵀmat(qb)

    return [X Q]
end
@inline function ∂g∂ʳposb(joint::AbstractJoint, xb::AbstractVector, qb::UnitQuaternion)
    X, Q = ∂g∂posb(joint, xb, qb)
    Q = Q * LVᵀmat(qb)

    return [X Q]
end

## Discrete-time velocity derivatives (for dynamics)
# Wrappers 1
@inline function ∂g∂ʳvela(joint::AbstractJoint, body1::Body, body2::Body, childid, Δt)
    if body2.id == childid
        return constraintmat(joint) * ∂g∂ʳvela(joint, body1.state, body2.state, Δt)
    else
        return zero(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::AbstractJoint, body1::Body, body2::Body, childid, Δt)
    if body2.id == childid
        return constraintmat(joint) * ∂g∂ʳvelb(joint, body1.state, body2.state, Δt)
    else
        return zero(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::AbstractJoint, body1::Origin, body2::Body, childid, Δt)
    if body2.id == childid
        return constraintmat(joint) * ∂g∂ʳvelb(joint, body2.state, Δt)
    else
        return zero(joint)
    end
end

@inline function offdiagonal∂damper∂ʳvel(joint::AbstractJoint, body1::Body, body2::Body, childid)
    if body2.id == childid
        return offdiagonal∂damper∂ʳvel(joint, body1.state, body2.state)
    else
        return zero(body2)
    end
end
@inline function offdiagonal∂damper∂ʳvel(joint::AbstractJoint, body1::Origin, body2::Body, childid)
    if body2.id == childid
        return offdiagonal∂damper∂ʳvel(joint, body2.state)
    else
        return zero(body2)
    end
end




# Wrappers 2
∂g∂ʳvela(joint::AbstractJoint, statea::State, stateb::State, Δt) = ∂g∂ʳvela(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)..., fullargssol(statea)..., Δt)
∂g∂ʳvelb(joint::AbstractJoint, statea::State, stateb::State, Δt) = ∂g∂ʳvelb(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)..., fullargssol(stateb)..., Δt)
∂g∂ʳvelb(joint::AbstractJoint, stateb::State, Δt) = ∂g∂ʳvelb(joint, posargsnext(stateb, Δt)..., fullargssol(stateb)..., Δt)
offdiagonal∂damper∂ʳvel(joint::AbstractJoint, statea::State, stateb::State) = offdiagonal∂damper∂ʳvel(joint, posargsk(statea)..., posargsk(stateb)...)
offdiagonal∂damper∂ʳvel(joint::AbstractJoint, stateb::State) = offdiagonal∂damper∂ʳvel(joint, posargsk(stateb)...)

# Derivatives accounting for quaternion specialness
@inline function ∂g∂ʳvela(joint::AbstractJoint, x2a::AbstractVector, q2a::UnitQuaternion, x2b::AbstractVector, q2b::UnitQuaternion,
        x1a::AbstractVector, v1a::AbstractVector, q1a::UnitQuaternion, ω1a::AbstractVector, Δt
    )

    X, Q = ∂g∂posa(joint, x2a, q2a, x2b, q2b)
    V = X * Δt
    Ω = Q * Lmat(q1a) * derivωbar(ω1a, Δt) * Δt / 2

    return [V Ω]
end
@inline function ∂g∂ʳvelb(joint::AbstractJoint, x2a::AbstractVector, q2a::UnitQuaternion, x2b::AbstractVector, q2b::UnitQuaternion,
        x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, Δt
    )

    X, Q = ∂g∂posb(joint, x2a, q2a, x2b, q2b)
    V = X * Δt
    Ω = Q * Lmat(q1b) * derivωbar(ω1b, Δt) * Δt / 2

    return [V Ω]
end
@inline function ∂g∂ʳvelb(joint::AbstractJoint, x2b::AbstractVector, q2b::UnitQuaternion,
        x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, Δt
    )

    X, Q = ∂g∂posb(joint, x2b, q2b)
    V = X * Δt
    Ω = Q * Lmat(q1b) * derivωbar(ω1b, Δt) * Δt / 2

    return [V Ω]
end

## Continuous-time position derivatives (for constraint solver)
# Wrappers 1
@inline function ∂g∂posac(joint::AbstractJoint, body1::Body, body2::Body, childid)
    if body2.id == childid
        return constraintmat(joint) * ∂g∂posac(joint, body1.state, body2.state)
    else
        return ∂g∂posac(joint)
    end
end
@inline function ∂g∂posbc(joint::AbstractJoint, body1::Body, body2::Body, childid)
    if body2.id == childid
        return constraintmat(joint) * ∂g∂posbc(joint, body1.state, body2.state)
    else
        return ∂g∂posbc(joint)
    end
end
@inline function ∂g∂posbc(joint::AbstractJoint, body1::Origin, body2::Body, childid)
    if body2.id == childid
        return constraintmat(joint) * ∂g∂posbc(joint, body2.state)
    else
        return ∂g∂posbc(joint)
    end
end

# Wrappers 2
@inline ∂g∂posac(joint::AbstractJoint{T,N}) where {T,N} = szeros(T, N, 7) # TODO zero function?
@inline ∂g∂posbc(joint::AbstractJoint{T,N}) where {T,N} = szeros(T, N, 7)
∂g∂posac(joint::AbstractJoint, statea::State, stateb::State) = hcat(∂g∂posa(joint, posargsc(statea)..., posargsc(stateb)...)...)
∂g∂posbc(joint::AbstractJoint, statea::State, stateb::State) = hcat(∂g∂posb(joint, posargsc(statea)..., posargsc(stateb)...)...)
∂g∂posbc(joint::AbstractJoint, stateb::State) = hcat(∂g∂posb(joint, posargsc(stateb)...)...)


### Springs and Dampers (for dynamics)
## Wrappers 1
@inline function springforcea(joint::AbstractJoint, body1::Body, body2::Body, childid)
    if body2.id == childid
        return springforcea(joint, body1.state, body2.state)
    else
        return springforce(joint)
    end
end
@inline function springforceb(joint::AbstractJoint, body1::Body, body2::Body, childid)
    if body2.id == childid
        return springforceb(joint, body1.state, body2.state)
    else
        return springforce(joint)
    end
end
@inline function springforceb(joint::AbstractJoint, body1::Origin, body2::Body, childid)
    if body2.id == childid
        return springforceb(joint, body2.state)
    else
        return springforce(joint)
    end
end

@inline function damperforcea(joint::AbstractJoint, body1::Body, body2::Body, childid)
    if body2.id == childid
        return damperforcea(joint, body1.state, body2.state)
    else
        return damperforce(joint)
    end
end
@inline function damperforceb(joint::AbstractJoint, body1::Body, body2::Body, childid)
    if body2.id == childid
        return damperforceb(joint, body1.state, body2.state)
    else
        return damperforce(joint)
    end
end
@inline function damperforceb(joint::AbstractJoint, body1::Origin, body2::Body, childid)
    if body2.id == childid
        return damperforceb(joint, body2.state)
    else
        return damperforce(joint)
    end
end

## Wrappers 2
@inline springforce(joint::AbstractJoint{T}) where {T} = szeros(T, 6) # TODO zero function?
springforcea(joint::AbstractJoint, statea::State, stateb::State) = springforcea(joint, posargsk(statea)..., posargsk(stateb)...)
springforceb(joint::AbstractJoint, statea::State, stateb::State) = springforceb(joint, posargsk(statea)..., posargsk(stateb)...)
springforceb(joint::AbstractJoint, stateb::State) = springforceb(joint, posargsk(stateb)...)

@inline damperforce(joint::AbstractJoint{T}) where {T} = szeros(T, 6) # TODO zero function?
damperforcea(joint::AbstractJoint, statea::State, stateb::State) = damperforcea(joint, fullargssol(statea)..., fullargssol(stateb)...)
damperforceb(joint::AbstractJoint, statea::State, stateb::State) = damperforceb(joint, fullargssol(statea)..., fullargssol(stateb)...)
damperforceb(joint::AbstractJoint, stateb::State) = damperforceb(joint, fullargssol(stateb)...)


### Forcing (for dynamics)
## Wrappers
@inline function applyFτ!(joint::AbstractJoint, body1::Body, body2::Body, clear::Bool)
    applyFτ!(joint, body1.state, body2.state, clear)
    return
end
@inline function applyFτ!(joint::AbstractJoint, ::Origin, body2::Body, clear::Bool)
    applyFτ!(joint, body2.state, clear)
    return
end
