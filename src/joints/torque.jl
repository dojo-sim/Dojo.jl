mutable struct Torque{T,N} <: FJoint{T,N}
    V3::Adjoint{T,SVector{3,T}} # in body1's frame
    V12::SMatrix{2,3,T,6} # in body1's frame
    vertices::NTuple{2,SVector{3,T}} # in body1's & body2's frames
    qoffset::UnitQuaternion{T} # in body1's frame

    spring::T
    damper::T

    #TODO: reference

    Fτ::SVector{3,T}

    function Torque{T,N}(body1::AbstractBody, body2::AbstractBody;
            p1::AbstractVector = szeros(T,3), p2::AbstractVector = szeros(T,3), axis::AbstractVector = szeros(T,3), qoffset::UnitQuaternion = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)
        ) where {T,N}

        vertices = (p1, p2)
        V1, V2, V3 = orthogonalrows(axis)
        V12 = [V1;V2]

        Fτ = zeros(T,3)

        new{T,N}(V3, V12, vertices, qoffset, spring, damper, Fτ), body1.id, body2.id
    end
end

Torque0 = Torque{T,0} where T
Torque1 = Torque{T,1} where T
Torque2 = Torque{T,2} where T
Torque3 = Torque{T,3} where T

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::Torque{T,N}) where {T,N}
    summary(io, constraint)
    println(io,"")
    println(io, " V3:       "*string(constraint.V3))
    println(io, " V12:      "*string(constraint.V12))
    println(io, " vertices: "*string(constraint.vertices))
end

### Constraints and derivatives
## Position level constraint wrappers
@inline function g(joint::Torque, body1::Body, body2::Body, Δt)
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    Fτ = springtorque(joint, body1.state, body2.state, Δt) + dampertorque(joint, body1.state, body2.state)
    return Aᵀ * A * Fτ
end

@inline function g(joint::Torque, body1::Origin, body2::Body, Δt)
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    Fτ = springtorque(joint, body2.state, Δt) + dampertorque(joint, body2.state)
    return Aᵀ * A * Fτ
end


## Position level constraints (for dynamics) in world frame
@inline function gc(joint::Torque, qa::UnitQuaternion, qb::UnitQuaternion)
    return Vmat(qa \ qb / joint.qoffset)
end
@inline function gc(joint::Torque, qb::UnitQuaternion)
    return Vmat(qb / joint.qoffset)
end


### Spring and damper
## Discrete-time position wrappers (for dynamics)
springtorque(joint::Torque, statea::State, stateb::State, Δt) = springtorque(joint, posargsk(statea)[2], posargsnext(statea, Δt)[2], posargsnext(stateb, Δt)[2])
springtorque(joint::Torque, stateb::State, Δt) = springtorque(joint, posargsnext(stateb, Δt)[2])

dampertorque(joint::Torque, statea::State, stateb::State) = dampertorque(joint, posargsk(statea)[2], statea.ωsol[2], posargsk(stateb)[2], stateb.ωsol[2])
dampertorque(joint::Torque, stateb::State) = dampertorque(joint, posargsk(stateb)[2], stateb.ωsol[2])

### Spring and damper
# Force applied by body a on body b expressed in world frame
@inline function springtorque(joint::Torque, q1a::UnitQuaternion, qa::UnitQuaternion, qb::UnitQuaternion)
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    distance = A * gc(joint, qa, qb)
    force = -Aᵀ * A * joint.spring * Aᵀ * distance # force in offset frame

    force = vrotate(force, q1a * joint.qoffset) # rotate back to world frame
    return force
end
# Force applied by origin on body b expressed in world frame
@inline function springtorque(joint::Torque, qb::UnitQuaternion)
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    distance = A * gc(joint, qb)
    force = -Aᵀ * A * joint.spring * Aᵀ * distance
    force = vrotate(force, joint.qoffset) # rotate back to world frame
    return force
end
# Force applied by body a on body b expressed in world frame
@inline function dampertorque(joint::Torque, q1a::UnitQuaternion, ωa::AbstractVector, q1b::UnitQuaternion, ωb::AbstractVector)
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    qoffset = joint.qoffset
    velocity = A * (vrotate(ωb, q1a \ q1b / qoffset) - vrotate(ωa, inv(qoffset))) # in offset frame
    force = -2 * Aᵀ * A * joint.damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
    force = vrotate(force, q1a * qoffset) # rotate back to world frame
    return force
end
# Force applied by origin on body b expressed in world frame
@inline function dampertorque(joint::Torque, q1b::UnitQuaternion, ωb::AbstractVector)
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    qoffset = joint.qoffset
    velocity = A * vrotate(ωb, q1b / qoffset) # in offset frame
    force = -2 * Aᵀ * A * joint.damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
    force = vrotate(force, qoffset) # rotate back to world frame
    return force
end

# Wrappers 2
∂g∂ʳposa(joint::Torque, statea::State, stateb::State) = ∂g∂ʳposa(joint, posargsk(statea)..., posargsk(stateb)...)
∂g∂ʳposb(joint::Torque, statea::State, stateb::State) = ∂g∂ʳposb(joint, posargsk(statea)..., posargsk(stateb)...)
∂g∂ʳposb(joint::Torque, stateb::State) = ∂g∂ʳposb(joint, posargsk(stateb)...)

# Derivatives accounting for quaternion specialness
@inline function ∂g∂ʳposa(joint::Torque{T,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    rot = VRᵀmat(inv(qa)) * LVᵀmat(inv(qa))
    X = szeros(T, 3, 3) # accounts for the fact that λsol[2] holds the force applied by body a on body b.
    Q = - Aᵀ * A * rot
    return [X Q]
end
@inline function ∂g∂ʳposb(joint::Torque{T,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    rot = VRᵀmat(inv(qb)) * LVᵀmat(inv(qb))
    X = szeros(T, 3, 3) # accounts for the fact that λsol[2] holds the force applied by body a on body b.
    Q = Aᵀ * A * rot
    return [X Q]
end
@inline function ∂g∂ʳposb(joint::Torque{T,N}, xb::AbstractVector, qb::UnitQuaternion) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    rot = VRᵀmat(inv(qb)) * LVᵀmat(inv(qb))
    X = szeros(T, 3, 3) # accounts for the fact that λsol[2] holds the force applied by body a on body b.
    Q = Aᵀ * A * rot
    return [X Q]
end


## Derivatives NOT accounting for quaternion specialness
# THIS IS USED IN DATAMAT, IT HAS TO BE THE DERIVATIVE OF g WRT THE POS VARIABLES (X3, Q3)
@inline function ∂g∂posa(joint::Torque, body1::Body, body2::Body, Δt)
    X, Q = ∂g∂posa(posargsk(body1.state)[2], posargsnext(body1.state, Δt)[2], body1.state.ωsol[2], posargsk(body2.state)[2], posargsnext(body2.state, Δt)[2], body2.state.ωsol[2], Δt) # the Δt factor comes from g(joint::FJoint
    return Δt * X, Δt * Q
end
@inline function ∂g∂posa1(joint::Torque, body1::Body, body2::Body, Δt)
    X, Q = ∂g∂posa1(posargsk(body1.state)[2], posargsnext(body1.state, Δt)[2], body1.state.ωsol[2], posargsk(body2.state)[2], posargsnext(body2.state, Δt)[2], body2.state.ωsol[2], Δt) # the Δt factor comes from g(joint::FJoint
    return Δt * X, Δt * Q
end
@inline function ∂g∂posb(joint::Torque, body1::Body, body2::Body, Δt)
    X, Q = ∂g∂posb(joint, posargsk(body1.state)[2], posargsnext(body1.state, Δt)[2], body1.state.ωsol[2], posargsk(body2.state)[2], posargsnext(body2.state, Δt)[2], body2.state.ωsol[2], Δt) # the Δt factor comes from g(joint::FJoint
    return Δt * X, Δt * Q
end
@inline function ∂g∂posb1(joint::Torque, body1::Body, body2::Body, Δt)
    X, Q = ∂g∂posb1(joint, posargsk(body1.state)[2], posargsnext(body1.state, Δt)[2], body1.state.ωsol[2], posargsk(body2.state)[2], posargsnext(body2.state, Δt)[2], body2.state.ωsol[2], Δt) # the Δt factor comes from g(joint::FJoint
    return Δt * X, Δt * Q
end
@inline function ∂g∂posb(joint::Torque, body1::Origin, body2::Body, Δt)
    X, Q = ∂g∂posb(joint, posargsk(body2.state)[2], posargsnext(body2.state, Δt)[2], body2.state.ωsol[2], Δt) # the Δt factor comes from g(joint::FJoint
    return Δt * X, Δt * Q
end
@inline function ∂g∂posb1(joint::Torque, body1::Origin, body2::Body, Δt)
    X, Q = ∂g∂posb1(joint, posargsk(body2.state)[2], posargsnext(body2.state, Δt)[2], body2.state.ωsol[2], Δt) # the Δt factor comes from g(joint::FJoint
    return Δt * X, Δt * Q
end
@inline function ∂g∂posa(joint::Torque{T,N}, q1a::UnitQuaternion, qa::UnitQuaternion, ωa::AbstractVector, q1b::UnitQuaternion, qb::UnitQuaternion, ωb::AbstractVector, Δt) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    τ_spring = springtorque(joint, q1a, qa, qb)
    τ_damp = dampertorque(joint, q1a, ωa, q1b, ωb)
    qoffset = joint.qoffset
    Xdamp = szeros(T, 3, 3)
    Qdamp = szeros(T, 3, 3)
    Xspring = szeros(T, 3, 3)
    Qspring = -∂vrotate∂p(τ_spring, q1a * joint.qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * Rmat(qb * inv(qoffset)) * Tmat()
    X = Xdamp + Xspring
    Q = Qdamp + Qspring

    return X, Q
end
@inline function ∂g∂posa1(joint::Torque{T,N}, q1a::UnitQuaternion, qa::UnitQuaternion, ωa::AbstractVector, q1b::UnitQuaternion, qb::UnitQuaternion, ωb::AbstractVector, Δt) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    τ_spring = springtorque(joint, q1a, qa, qb)
    τ_damp = dampertorque(joint, q1a, ωa, q1b, ωb)
    qoffset = joint.qoffset

    Xdamp = szeros(T, 3, 3)
    Qdamp = ∂vrotate∂q(τ_damp, q1a * qoffset) * Rmat(qoffset)
    Qdamp += -2.0 * ∂vrotate∂p(τ_damp, q1a * qoffset) * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, q1a \ q1b / qoffset) * Rmat(q1b * inv(qoffset)) * Tmat()
    Xspring = szeros(T, 3, 3)
    Qspring = ∂vrotate∂q(τ_spring, q1a * qoffset) * Rmat(qoffset)
    X = Xdamp + Xspring
    Q = Qdamp + Qspring

    return X, Q
end
@inline function ∂g∂posb(joint::Torque{T,N}, q1a::UnitQuaternion, qa::UnitQuaternion, ωa::AbstractVector, q1b::UnitQuaternion, qb::UnitQuaternion, ωb::AbstractVector, Δt) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    τ_spring = springtorque(joint, q1a, qa, qb)
    τ_damp = dampertorque(joint, q1a, ωa, q1b, ωb)
    qoffset = joint.qoffset

    Xdamp = szeros(T, 3, 3)
    Qdamp = szeros(T, 3, 3)
    Xspring = szeros(T, 3, 3)
    Qspring = -1.0 * ∂vrotate∂p(τ_spring, q1a * qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * Lmat(inv(qa)) * Rmat(inv(qoffset))
    X = Xdamp + Xspring
    Q = Qdamp + Qspring

    return X, Q
end
@inline function ∂g∂posb1(joint::Torque{T,N}, q1a::UnitQuaternion, qa::UnitQuaternion, ωa::AbstractVector, q1b::UnitQuaternion, qb::UnitQuaternion, ωb::AbstractVector, Δt) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    τ_spring = springtorque(joint, q1a, qa, qb)
    τ_damp = dampertorque(joint, q1a, ωa, q1b, ωb)
    qoffset = joint.qoffset

    Xdamp = szeros(T, 3, 3)
    Qdamp = ∂vrotate∂p(τ_damp, q1a * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, q1a \ q1b / qoffset) * Lmat(inv(q1a)) * Rmat(inv(qoffset))
    Xspring = szeros(T, 3, 3)
    Qspring = szeros(T, 3, 3)
    X = Xdamp + Xspring
    Q = Qdamp + Qspring

    return X, Q
end
@inline function ∂g∂posb(joint::Torque{T,N}, q1b::UnitQuaternion, qb::UnitQuaternion, ωb::AbstractVector, Δt) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    τ_spring = springtorque(joint, qb)
    τ_damp = dampertorque(joint, q1b, ωb)
    qoffset = joint.qoffset

    Xdamp = szeros(T, 3, 3)
    Qdamp = szeros(T, 3, 4)
    Xspring = szeros(T, 3, 3)
    Qspring = -1.0 * ∂vrotate∂p(τ_spring, qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset))
    X = Xdamp + Xspring
    Q = Qdamp + Qspring

    return X, Q
end
@inline function ∂g∂posb1(joint::Torque{T,N}, q1b::UnitQuaternion, qb::UnitQuaternion, ωb::AbstractVector, Δt) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    τ_spring = springtorque(joint, qb)
    τ_damp = dampertorque(joint, q1b, ωb)
    qoffset = joint.qoffset

    Xdamp = szeros(T, 3, 3)
    Qdamp = -2 * ∂vrotate∂p(τ_damp, qoffset) * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, q1b / qoffset) * Rmat(inv(qoffset))
    Xspring = szeros(T, 3, 3)
    Qspring = szeros(T, 3, 4)
    X = Xdamp + Xspring
    Q = Qdamp + Qspring
    return X, Q
end

# Wrappers 2
∂g∂ʳvela(joint::Torque, statea::State, stateb::State, Δt) = ∂g∂ʳvela(joint, posargsc(statea)[2], posargsnext(statea, Δt)..., statea.vsol[2], statea.ωsol[2], posargsc(stateb)[2], posargsnext(stateb, Δt)..., stateb.vsol[2], stateb.ωsol[2], Δt)
∂g∂ʳvelb(joint::Torque, statea::State, stateb::State, Δt) = ∂g∂ʳvelb(joint, posargsc(statea)[2], posargsnext(statea, Δt)..., statea.vsol[2], statea.ωsol[2], posargsc(stateb)[2], posargsnext(stateb, Δt)..., stateb.vsol[2], stateb.ωsol[2], Δt)
∂g∂ʳvelb(joint::Torque, stateb::State, Δt) = ∂g∂ʳvelb(joint, posargsc(stateb)[2], posargsnext(stateb, Δt)..., stateb.vsol[2], stateb.ωsol[2], Δt)

# Derivatives accounting for quaternion specialness
@inline function ∂g∂ʳvela(joint::Torque{T,N}, q1a::UnitQuaternion, xa::AbstractVector,
        qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
        q1b::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion,
        vb::AbstractVector, ωb::AbstractVector, Δt) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    τ_spring = springtorque(joint, q1a, qa, qb)
    τ_damp = dampertorque(joint, q1a, ωa, q1b, ωb)
    qoffset = joint.qoffset
    Vdamp = szeros(T, 3, 3)
    Ωdamp = 2.0 * ∂vrotate∂p(τ_damp, q1a * qoffset) * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωa, inv(qoffset)) * Lmat(q1a) * derivωbar(ωa, Δt) * Δt/2
    Vspring = szeros(T, 3, 3)

    Ωspring = -1.0 * ∂vrotate∂p(τ_spring, q1a * joint.qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(qb * inv(qoffset)) * Tmat() * Lmat(q1a) * derivωbar(ωa, Δt) * Δt/2

    V = Vspring + Vdamp
    Ω = Ωspring + Ωdamp
    V *= Δt
    Ω *= Δt

    return [V Ω]
end
@inline function ∂g∂ʳvelb(joint::Torque{T,N}, q1a::UnitQuaternion, xa::AbstractVector,
        qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
        q1b::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion,
        vb::AbstractVector, ωb::AbstractVector, Δt) where {T,N}

    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    τ_spring = springtorque(joint, q1a, qa, qb)
    τ_damp = dampertorque(joint, q1a, ωa, q1b, ωb)
    qoffset = joint.qoffset
    Vdamp = szeros(T, 3, 3)
    Ωdamp = -2.0 * ∂vrotate∂p(τ_damp, q1a * qoffset) * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂p(ωb, q1a \ q1b / qoffset)
    Vspring = szeros(T, 3, 3)
    Ωspring = -1.0 * ∂vrotate∂p(τ_spring, q1a * qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VLmat(inv(qa)) * Rmat(inv(qoffset)) * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2

    V = Vspring + Vdamp
    Ω = Ωspring + Ωdamp
    V *= Δt
    Ω *= Δt

    return [V Ω]
end
@inline function ∂g∂ʳvelb(joint::Torque{T,N}, q1b::UnitQuaternion, xb::AbstractVector,
        qb::UnitQuaternion, vb::AbstractVector,  ωb::AbstractVector, Δt) where {T,N}
        A = constraintmat(joint)
        Aᵀ = zerodimstaticadjoint(A)
        τ_spring = springtorque(joint, qb)
        τ_damp = dampertorque(joint, q1b, ωb)
        qoffset = joint.qoffset
        Vdamp = szeros(T, 3, 3)
        Ωdamp = -2.0 * ∂vrotate∂p(τ_damp, qoffset) * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂p(ωb, q1b / qoffset)
        Vspring = szeros(T, 3, 3)
        Ωspring = -1.0 * ∂vrotate∂p(τ_spring, qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2

        V = Vspring + Vdamp
        Ω = Ωspring + Ωdamp
        V *= Δt
        Ω *= Δt
    return [V Ω]
end

## vec(G) Jacobian (also NOT accounting for quaternion specialness in the second derivative: ∂(∂ʳg∂posx)∂y)
@inline function ∂2g∂posaa(joint::Torque{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    Lpos = Lmat(UnitQuaternion(xb + vrotate(joint.vertices[2], qb) - xa))
    Ltpos = Lᵀmat(UnitQuaternion(xb + vrotate(joint.vertices[2], qb) - xa))

    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = szeros(T, 9, 4)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posab(joint::Torque{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = szeros(T, 9, 4)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posba(joint::Torque{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = szeros(T, 9, 4)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posbb(joint::Torque{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = szeros(T, 9, 4)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posbb(joint::Torque{T}, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = szeros(T, 9, 4)

    return XX, XQ, QX, QQ
end

### Forcing
## Application of joint forces (for dynamics)
@inline function applyFτ!(joint::Torque{T}, statea::State, stateb::State, clear::Bool) where T
    τ = joint.Fτ
    _, qa = posargsk(statea)
    _, qb = posargsk(stateb)

    τa = -τ # in world coordinates
    τb = τ  # in world coordinates

    τa = vrotate(τa, inv(qa)) # in local coordinates
    τb = vrotate(τb, inv(qb)) # in local coordinates

    statea.τk[end] += τa
    stateb.τk[end] += τb
    clear && (joint.Fτ = szeros(T,3))
    return
end
@inline function applyFτ!(joint::Torque{T}, stateb::State, clear::Bool) where T
    τ = joint.Fτ
    _, qb = posargsk(stateb)

    τb = τ # in world coordinates

    τb = vrotate(τb,inv(qb)) # in local coordinates

    stateb.τk[end] += τb
    clear && (joint.Fτ = szeros(T,3))
    return
end

## Forcing derivatives (for linearization)
# Control derivatives
@inline function ∂Fτ∂ua(joint::Torque, statea::State, stateb::State)
    BFa = szeros(T, 3, 3)
    Bτa = szeros(T, 3, 3)

    return [BFa; Bτa]
end

@inline function ∂Fτ∂ub(joint::Torque, statea::State, stateb::State)
    BFb = szeros(T, 3, 3)
    Bτb = szeros(T, 3, 3)

    return [BFb; Bτb]
end
@inline function ∂Fτ∂ub(joint::Torque, stateb::State)
    BFb = szeros(T, 3, 3)
    Bτb = szeros(T, 3, 3)
    return [BFb; Bτb]
end

# Position derivatives
@inline function ∂Fτ∂posa(joint::Torque{T}, statea::State, stateb::State) where T
    FaXa = szeros(T,3,3)
    FaQa = szeros(T,3,3)
    τaXa = szeros(T,3,3)
    τaQa = szeros(T,3,3)
    FbXa = szeros(T,3,3)
    FbQa = szeros(T,3,3)
    τbXa = szeros(T,3,3)
    τbQa = szeros(T,3,3)

    return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
end
@inline function ∂Fτ∂posb(joint::Torque{T}, statea::State, stateb::State) where T
    FaXb = szeros(T,3,3)
    FaQb = szeros(T,3,3)
    τaXb = szeros(T,3,3)
    τaQb = szeros(T,3,3)
    FbXb = szeros(T,3,3)
    FbQb = szeros(T,3,3)
    τbXb = szeros(T,3,3)
    τbQb = szeros(T,3,3)

    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end
@inline function ∂Fτ∂posb(joint::Torque{T}, stateb::State) where T
    xb, qb = posargsk(stateb)
    F = joint.Fτ
    vertices = joint.vertices

    FaXb = szeros(T,3,3)
    FaQb = szeros(T,3,3)
    τaXb = szeros(T,3,3)
    τaQb = szeros(T,3,3)
    FbXb = szeros(T,3,3)
    FbQb = szeros(T,3,3)
    τbXb = szeros(T,3,3)
    τbQb = szeros(T,3,3)

    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end


# ### Minimal coordinates
# ## Position and velocity offsets
# @inline function getPositionDelta(joint::Torque, body1::AbstractBody, body2::Body, x::SVector)
#     Δx = zerodimstaticadjoint(nullspacemat(joint)) * x # in body1 frame
#     return Δx
# end
# @inline function getVelocityDelta(joint::Torque, body1::AbstractBody, body2::Body, v::SVector)
#     Δv = zerodimstaticadjoint(nullspacemat(joint)) * v # in body1 frame
#     return Δv
# end

# ## Minimal coordinate calculation
# @inline function minimalCoordinates(joint::Torque, body1::Body, body2::Body)
#     statea = body1.state
#     stateb = body2.state
#     return nullspacemat(joint) * g(joint, statea.xc, statea.qc, stateb.xc, stateb.qc)
# end
# @inline function minimalCoordinates(joint::Torque, body1::Origin, body2::Body)
#     stateb = body2.state
#     return nullspacemat(joint) * g(joint, stateb.xc, stateb.qc)
# end
# @inline function minimalVelocities(joint::Torque, body1::Body, body2::Body)
#     statea = body1.state
#     stateb = body2.state
#     return nullspacemat(joint) * (stateb.vc - statea.vc)
# end
# @inline function minimalVelocities(joint::Torque, body1::Origin, body2::Body)
#     stateb = body2.state
#     return nullspacemat(joint) * stateb.vc
# end
# nullspacemat(force1)
