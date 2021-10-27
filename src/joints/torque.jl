# mutable struct Rotational{T,N} <: FJoint{T,N}
#     V3::Adjoint{T,SVector{3,T}} # in body1's frame
#     V12::SMatrix{2,3,T,6} # in body1's frame
#     vertices::NTuple{2,SVector{3,T}} # in body1's & body2's frames
#     qoffset::UnitQuaternion{T} # in body1's frame
#
#     spring::T
#     damper::T
#
#     #TODO: reference
#
#     Fτ::SVector{3,T}
#
#     function Rotational{T,N}(body1::AbstractBody, body2::AbstractBody;
#             p1::AbstractVector = szeros(T,3), p2::AbstractVector = szeros(T,3), axis::AbstractVector = szeros(T,3), qoffset::UnitQuaternion = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)
#         ) where {T,N}
#
#         vertices = (p1, p2)
#         V1, V2, V3 = orthogonalrows(axis)
#         V12 = [V1;V2]
#
#         Fτ = zeros(T,3)
#
#         new{T,N}(V3, V12, vertices, qoffset, spring, damper, Fτ), body1.id, body2.id
#     end
# end
#
# Rotational0 = Rotational{T,0} where T
# Rotational1 = Rotational{T,1} where T
# Rotational2 = Rotational{T,2} where T
# Rotational3 = Rotational{T,3} where T
#
# function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::Rotational{T,N}) where {T,N}
#     summary(io, constraint)
#     println(io,"")
#     println(io, " V3:       "*string(constraint.V3))
#     println(io, " V12:      "*string(constraint.V12))
#     println(io, " vertices: "*string(constraint.vertices))
# end

# ### Constraints and derivatives
# ## Position level constraint wrappers
# @inline function g(joint::Rotational, body1::Body, body2::Body, Δt)
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#     Fτ = springtorque(joint, body1.state, body2.state, Δt) + dampertorque(joint, body1.state, body2.state, Δt)
#     return Aᵀ * A * Fτ
# end
#
# @inline function g(joint::Rotational, body1::Origin, body2::Body, Δt)
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#     Fτ = springtorque(joint, body2.state, Δt) + dampertorque(joint, body2.state, Δt)
#     return Aᵀ * A * Fτ
# end

#
## Position level constraints (for dynamics) in world frame
@inline function gc(joint::Rotational, qa::UnitQuaternion, qb::UnitQuaternion)
    return Vmat(qa \ qb / joint.qoffset)
end
@inline function gc(joint::Rotational, qb::UnitQuaternion)
    return Vmat(qb / joint.qoffset)
end


### Spring and damper
## Discrete-time position wrappers (for dynamics)
springtorque(joint::Rotational, statea::State, stateb::State, Δt) = springtorque(joint, posargsnext(statea, Δt)[2], posargsnext(stateb, Δt)[2])
springtorque(joint::Rotational, stateb::State, Δt) = springtorque(joint, posargsnext(stateb, Δt)[2])

dampertorque(joint::Rotational, statea::State, stateb::State, Δt) = dampertorque(joint, posargsnext(statea, Δt)[2], statea.ωsol[2], posargsnext(stateb, Δt)[2], stateb.ωsol[2])
dampertorque(joint::Rotational, stateb::State, Δt) = dampertorque(joint, posargsnext(stateb, Δt)[2], stateb.ωsol[2])

### Spring and damper
# Force applied by body b on body a expressed in frame a
@inline function springtorquea(joint::Rotational, qa::UnitQuaternion, qb::UnitQuaternion; rotate::Bool = true)
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    distance = A * gc(joint, qa, qb)
    force = Aᵀ * A * joint.spring * Aᵀ * distance # force in offset frame
    rotate && (force = vrotate(force, joint.qoffset)) # rotate back to a frame
    return force
end
# Force applied by body a on body b expressed in frame b
@inline function springtorqueb(joint::Rotational, qa::UnitQuaternion, qb::UnitQuaternion; rotate::Bool = true)
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    distance = A * gc(joint, qa, qb)
    force = - Aᵀ * A * joint.spring * Aᵀ * distance # force in offset frame
    rotate && (force = vrotate(force, inv(qb) * qa * joint.qoffset)) # rotate back to b frame
    return force
end
# Force applied by origin on body b expressed in frame b
@inline function springtorqueb(joint::Rotational, qb::UnitQuaternion; rotate::Bool = true)
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    distance = A * gc(joint, qb)
    force = - Aᵀ * A * joint.spring * Aᵀ * distance # force in offset frame
    rotate && (force = vrotate(force, inv(qb) * joint.qoffset)) # rotate back to b frame
    return force
end

# Force applied by body b on body a expressed in frame a
@inline function dampertorquea(joint::Rotational, qa::UnitQuaternion, ωa::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector; rotate::Bool = true)
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    qoffset = joint.qoffset
    velocity = A * (vrotate(ωb, qa \ qb / qoffset) - vrotate(ωa, inv(qoffset))) # in offset frame
    force = 2 * Aᵀ * A * joint.damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
    rotate && (force = vrotate(force, qoffset)) # rotate back to frame a
    return force
end
# Force applied by body a on body b expressed in frame b
@inline function dampertorqueb(joint::Rotational, qa::UnitQuaternion, ωa::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector; rotate::Bool = true)
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    qoffset = joint.qoffset
    velocity = A * (vrotate(ωb, qa \ qb / qoffset) - vrotate(ωa, inv(qoffset))) # in offset frame
    force = - 2 * Aᵀ * A * joint.damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
    rotate && (force = vrotate(force, inv(qb) * qa * qoffset)) # rotate back to frame b
    return force
end
# Force applied by origin on body b expressed in frame b
@inline function dampertorqueb(joint::Rotational, qb::UnitQuaternion, ωb::AbstractVector; rotate::Bool = true)
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    qoffset = joint.qoffset
    velocity = A * vrotate(ωb, qb / qoffset) # in offset frame
    force = - 2 * Aᵀ * A * joint.damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
    rotate && (force = vrotate(force, inv(qb) * qoffset)) # rotate back to frame b
    return force
end

function ∂springtorquea∂posa(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1a, q1a = posargsk(body1.state)
    _, _, _, ωa = fullargssol(body1.state)
    xa, qa = posargsnext(body1.state, Δt)
    xb, qb = posargsnext(body2.state, Δt)
    qoffset = joint.qoffset
    force = springtorquea(joint, qa, qb; rotate = false)
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(qb * inv(qoffset)) * Tmat() * Rmat(ωbar(ωa, Δt)*Δt/2) * LVᵀmat(q1a)
    return [szeros(T, 3, 6); X Q]
end
function ∂dampertorquea∂posa(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1a, q1a = posargsk(body1.state)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargsnext(body1.state, Δt)
    xb, qb = posargsnext(body2.state, Δt)
    qoffset = joint.qoffset
    force = dampertorquea(joint, qa, ωa, qb, ωb; rotate = false)
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, qoffset) * 2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(qb * inv(qoffset)) * Tmat() * Rmat(ωbar(ωa, Δt)*Δt/2) * LVᵀmat(q1a)
    return [szeros(T, 3, 6); X Q]
end
function ∂springtorquea∂posb(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1b, q1b = posargsk(body2.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargsnext(body1.state, Δt)
    xb, qb = posargsnext(body2.state, Δt)
    qoffset = joint.qoffset
    force = springtorquea(joint, qa, qb; rotate = false)
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * Lmat(inv(qa)) * Rmat(ωbar(ωb, Δt)*Δt/2) * LVᵀmat(q1b)
    return [szeros(T, 3, 6); X Q]
end
function ∂dampertorquea∂posb(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1b, q1b = posargsk(body2.state)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargsnext(body1.state, Δt)
    xb, qb = posargsnext(body2.state, Δt)
    qoffset = joint.qoffset
    force = dampertorquea(joint, qa, ωa, qb, ωb; rotate = false)
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, qoffset) * 2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(inv(qoffset)) * Lmat(inv(qa)) * Rmat(ωbar(ωb, Δt)*Δt/2) * LVᵀmat(q1b)
    return [szeros(T, 3, 6); X Q]
end
function ∂springtorqueb∂posb(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1b, q1b = posargsk(body2.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargsnext(body1.state, Δt)
    xb, qb = posargsnext(body2.state, Δt)
    qoffset = joint.qoffset
    force = springtorqueb(joint, qa, qb; rotate = false)
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -1.0 * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * Lmat(inv(qa)) * Rmat(ωbar(ωb, Δt)*Δt/2) * LVᵀmat(q1b)
    Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qa * qoffset) * Tmat() * Rmat(ωbar(ωb, Δt)*Δt/2) * LVᵀmat(q1b)
    return [szeros(T, 3, 6); X Q]
end
function ∂dampertorqueb∂posb(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1b, q1b = posargsk(body2.state)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargsnext(body1.state, Δt)
    xb, qb = posargsnext(body2.state, Δt)
    qoffset = joint.qoffset
    force = dampertorqueb(joint, qa, ωa, qb, ωb; rotate = false)
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(inv(qoffset)) * Lmat(inv(qa)) * Rmat(ωbar(ωb, Δt)*Δt/2) * LVᵀmat(q1b)
    Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qa * qoffset) * Tmat() * Rmat(ωbar(ωb, Δt)*Δt/2) * LVᵀmat(q1b)
    return [szeros(T, 3, 6); X Q]
end
function ∂springtorqueb∂posa(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1a, q1a = posargsk(body1.state)
    _, _, _, ωa = fullargssol(body1.state)
    xa, qa = posargsnext(body1.state, Δt)
    xb, qb = posargsnext(body2.state, Δt)
    qoffset = joint.qoffset
    force = springtorqueb(joint, qa, qb; rotate = false)
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -1.0 * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(qb * inv(qoffset)) * Tmat() * Rmat(ωbar(ωa, Δt)*Δt/2) * LVᵀmat(q1a)
    Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qoffset) * Lmat(inv(qb)) * Rmat(ωbar(ωa, Δt)*Δt/2) * LVᵀmat(q1a)
    return [szeros(T, 3, 6); X Q]
end
function ∂dampertorqueb∂posa(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1a, q1a = posargsk(body1.state)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargsnext(body1.state, Δt)
    xb, qb = posargsnext(body2.state, Δt)
    qoffset = joint.qoffset
    force = dampertorqueb(joint, qa, ωa, qb, ωb; rotate = false)
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(qb * inv(qoffset)) * Tmat() * Rmat(ωbar(ωa, Δt)*Δt/2) * LVᵀmat(q1a)
    Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qoffset) * Lmat(inv(qb)) * Rmat(ωbar(ωa, Δt)*Δt/2) * LVᵀmat(q1a)
    return [szeros(T, 3, 6); X Q]
end
function ∂springtorqueb∂posb(joint::Rotational, body1::Origin, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1b, q1b = posargsk(body2.state)
    _, _, _, ωb = fullargssol(body2.state)
    xb, qb = posargsnext(body2.state, Δt)
    qoffset = joint.qoffset
    force = springtorqueb(joint, qb; rotate = false)
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, inv(qb) * qoffset) * -1.0 * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * Rmat(ωbar(ωb, Δt)*Δt/2) * LVᵀmat(q1b)
    Q += ∂vrotate∂q(force, inv(qb) * qoffset) * Rmat(qoffset) * Tmat() * Rmat(ωbar(ωb, Δt)*Δt/2) * LVᵀmat(q1b)
    return [szeros(T, 3, 6); X Q]
end
function ∂dampertorqueb∂posb(joint::Rotational, body1::Origin, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1b, q1b = posargsk(body2.state)
    _, _, _, ωb = fullargssol(body2.state)
    xb, qb = posargsnext(body2.state, Δt)
    qoffset = joint.qoffset
    force = dampertorqueb(joint, qb, ωb; rotate = false)
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, inv(qb) * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qb / qoffset) * Rmat(inv(qoffset)) * Rmat(ωbar(ωb, Δt)*Δt/2) * LVᵀmat(q1b)
    Q += ∂vrotate∂q(force, inv(qb) * qoffset) * Rmat(qoffset) * Tmat() * Rmat(ωbar(ωb, Δt)*Δt/2) * LVᵀmat(q1b)
    return [szeros(T, 3, 6); X Q]
end


function ∂springtorquea∂vela(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1a, q1a = posargsk(body1.state)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargsnext(body1.state, Δt)
    xb, qb = posargsnext(body2.state, Δt)
    qoffset = joint.qoffset
    force = springtorquea(joint, qa, qb; rotate = false)
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(qb * inv(qoffset)) * Tmat() * Lmat(q1a) * derivωbar(ωa, Δt) * Δt/2
    return [szeros(T, 3, 6); V Ω]
end
function ∂dampertorquea∂vela(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    V = - Aᵀ * A * joint.damper * Aᵀ * A
    Ω = szeros(T, 3, 3)
    return [szeros(T, 3, 6); V Ω]
end
function ∂springtorquea∂velb(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1b, q1b = posargsk(body2.state)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargsnext(body1.state, Δt)
    xb, qb = posargsnext(body2.state, Δt)
    qoffset = joint.qoffset
    force = springtorquea(joint, qa, qb; rotate = false)
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * Lmat(inv(qa)) * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2
    return [szeros(T, 3, 6); V Ω]
end
function ∂dampertorquea∂velb(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    V = Aᵀ * A * joint.damper * Aᵀ * A
    Ω = szeros(T, 3, 3)
    return [szeros(T, 3, 6); V Ω]
end
function ∂springtorqueb∂velb(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1b, q1b = posargsk(body2.state)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargsnext(body1.state, Δt)
    xb, qb = posargsnext(body2.state, Δt)
    qoffset = joint.qoffset
    force = springtorqueb(joint, qa, qb; rotate = false)
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -1.0 * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * Lmat(inv(qa)) * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2
    Ω += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qa * qoffset) * Tmat() * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2
    return [szeros(T, 3, 6); V Ω]
end
function ∂dampertorqueb∂velb(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    V = - Aᵀ * A * joint.damper * Aᵀ * A
    Ω = szeros(T, 3, 3)
    return [szeros(T, 3, 6); V Ω]
end
function ∂springtorqueb∂vela(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1a, q1a = posargsk(body1.state)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargsnext(body1.state, Δt)
    xb, qb = posargsnext(body2.state, Δt)
    qoffset = joint.qoffset
    force = springtorqueb(joint, qa, qb; rotate = false)
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -1.0 * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(qb * inv(qoffset)) * Tmat() * Lmat(q1a) * derivωbar(ωa, Δt) * Δt/2
    Ω += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qoffset) * Lmat(inv(qb)) * Lmat(q1a) * derivωbar(ωa, Δt) * Δt/2
    return [szeros(T, 3, 6); V Ω]
end
function ∂dampertorqueb∂vela(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    V = Aᵀ * A * joint.damper * Aᵀ * A
    Ω = szeros(T, 3, 3)
    return [szeros(T, 3, 6); V Ω]
end
function ∂springtorqueb∂velb(joint::Rotational, body1::Origin, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1b, q1b = posargsk(body2.state)
    _, _, _, ωb = fullargssol(body2.state)
    xb, qb = posargsnext(body2.state, Δt)
    qoffset = joint.qoffset
    force = springtorqueb(joint, qb; rotate = false)
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, inv(qb) * qoffset) * -1.0 * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2
    Ω += ∂vrotate∂q(force, inv(qb) * qoffset) * Rmat(qoffset) * Tmat() * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2
    return [szeros(T, 3, 6); V Ω]
end
function ∂dampertorqueb∂velb(joint::Rotational, body1::Origin, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    V = - Aᵀ * A * joint.damper * Aᵀ * A
    Ω = szeros(T, 3, 3)
    return [szeros(T, 3, 6); V Ω]
end



# # Wrappers 2
# ∂g∂ʳposa(joint::Rotational, statea::State, stateb::State, Δt) = ∂g∂ʳposa(joint, posargsk(statea)..., posargsk(stateb)...)
# ∂g∂ʳposb(joint::Rotational, statea::State, stateb::State, Δt) = ∂g∂ʳposb(joint, posargsk(statea)..., posargsk(stateb)...)
# ∂g∂ʳposb(joint::Rotational, stateb::State, Δt) = ∂g∂ʳposb(joint, posargsk(stateb)...)
#
# # Derivatives accounting for quaternion specialness
# @inline function ∂g∂ʳposa(joint::Rotational{T,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where {T,N}
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#     rot = VRᵀmat(inv(qa)) * LVᵀmat(inv(qa))
#     X = szeros(T, 3, 3) # accounts for the fact that λsol[2] holds the force applied by body a on body b.
#     Q = - Aᵀ * A * rot
#     return [X Q]
# end
# @inline function ∂g∂ʳposb(joint::Rotational{T,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where {T,N}
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#     rot = VRᵀmat(inv(qb)) * LVᵀmat(inv(qb))
#     X = szeros(T, 3, 3) # accounts for the fact that λsol[2] holds the force applied by body a on body b.
#     Q = Aᵀ * A * rot
#     return [X Q]
# end
# @inline function ∂g∂ʳposb(joint::Rotational{T,N}, xb::AbstractVector, qb::UnitQuaternion) where {T,N}
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#     rot = VRᵀmat(inv(qb)) * LVᵀmat(inv(qb))
#     X = szeros(T, 3, 3) # accounts for the fact that λsol[2] holds the force applied by body a on body b.
#     Q = Aᵀ * A * rot
#     return [X Q]
# end
#
#
# ## Derivatives NOT accounting for quaternion specialness
# # THIS IS USED IN DATAMAT, IT HAS TO BE THE DERIVATIVE OF g WRT THE POS VARIABLES (X3, Q3)
# @inline function ∂g∂posa(joint::Rotational, body1::Body, body2::Body, Δt)
#     X, Q = ∂g∂posa(joint, posargsnext(body1.state, Δt)[2], body1.state.ωsol[2], posargsnext(body2.state, Δt)[2], body2.state.ωsol[2], Δt) # the Δt factor comes from g(joint::FJoint
#     return Δt * X, Δt * Q
# end
# @inline function ∂g∂posb(joint::Rotational, body1::Body, body2::Body, Δt)
#     X, Q = ∂g∂posb(joint, posargsnext(body1.state, Δt)[2], body1.state.ωsol[2], posargsnext(body2.state, Δt)[2], body2.state.ωsol[2], Δt) # the Δt factor comes from g(joint::FJoint
#     return Δt * X, Δt * Q
# end
# @inline function ∂g∂posb(joint::Rotational, body1::Origin, body2::Body, Δt)
#     X, Q = ∂g∂posb(joint, posargsnext(body2.state, Δt)[2], body2.state.ωsol[2], Δt) # the Δt factor comes from g(joint::FJoint
#     return Δt * X, Δt * Q
# end
@inline function ∂g∂posa(joint::Rotational{T,N}, qa::UnitQuaternion, ωa::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector, Δt) where {T,N}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    qoffset = joint.qoffset
    τ_spring = springtorque(joint, qa, qb, rotate = false)
    τ_damp = dampertorque(joint, qa, ωa, qb, ωb, rotate = false)

    Xdamp = szeros(T, 3, 3)
    Qdamp = ∂vrotate∂p(τ_damp, qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(qb * inv(qoffset)) * Tmat()
    Qdamp += ∂vrotate∂q(τ_damp, qa * qoffset) * Rmat(qoffset)
    Xspring = szeros(T, 3, 3)
    Qspring = ∂vrotate∂p(τ_spring, qa * qoffset) * -Aᵀ * A * joint.spring * Aᵀ * A * VRmat(qb * inv(qoffset)) * Tmat()
    Qspring += ∂vrotate∂q(τ_spring, qa * qoffset) * Rmat(qoffset)
    X = Xdamp + Xspring
    Q = Qdamp + Qspring

    return Aᵀ * A * X, Aᵀ * A * Q
end
#
# @inline function ∂g∂posb(joint::Rotational{T,N}, qa::UnitQuaternion, ωa::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector, Δt) where {T,N}
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#     τ_spring = springtorque(joint, qa, qb, rotate = false)
#     τ_damp = dampertorque(joint, qa, ωa, qb, ωb, rotate = false)
#     qoffset = joint.qoffset
#
#     Xdamp = szeros(T, 3, 3)
#     Qdamp = ∂vrotate∂p(τ_damp, qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(inv(qoffset)) * Lmat(inv(qa))
#     Xspring = szeros(T, 3, 3)
#     Qspring = ∂vrotate∂p(τ_spring, qa * qoffset) * -Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * Lmat(inv(qa))
#     X = Xdamp + Xspring
#     Q = Qdamp + Qspring
#
#     return Aᵀ * A * X, Aᵀ * A * Q
# end
# @inline function ∂g∂posb(joint::Rotational{T,N}, qb::UnitQuaternion, ωb::AbstractVector, Δt) where {T,N}
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#     τ_spring = springtorque(joint, qb, rotate = false)
#     τ_damp = dampertorque(joint, qb, ωb, rotate = false)
#     qoffset = joint.qoffset
#
#     Xdamp = szeros(T, 3, 3)
#     Qdamp = ∂vrotate∂p(τ_damp, qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qb / qoffset) * Rmat(inv(qoffset))
#     Xspring = szeros(T, 3, 3)
#     Qspring = -1.0 * ∂vrotate∂p(τ_spring, qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset))
#     X = Xdamp + Xspring
#     Q = Qdamp + Qspring
#
#     return Aᵀ * A * X, Aᵀ * A * Q
# end
#
# # Wrappers 2
# ∂g∂ʳvela(joint::Rotational, statea::State, stateb::State, Δt) = ∂g∂ʳvela(joint, posargsc(statea)[2], posargsnext(statea, Δt)..., statea.vsol[2], statea.ωsol[2], posargsc(stateb)[2], posargsnext(stateb, Δt)..., stateb.vsol[2], stateb.ωsol[2], Δt)
# ∂g∂ʳvelb(joint::Rotational, statea::State, stateb::State, Δt) = ∂g∂ʳvelb(joint, posargsc(statea)[2], posargsnext(statea, Δt)..., statea.vsol[2], statea.ωsol[2], posargsc(stateb)[2], posargsnext(stateb, Δt)..., stateb.vsol[2], stateb.ωsol[2], Δt)
# ∂g∂ʳvelb(joint::Rotational, stateb::State, Δt) = ∂g∂ʳvelb(joint, posargsc(stateb)[2], posargsnext(stateb, Δt)..., stateb.vsol[2], stateb.ωsol[2], Δt)
# # Derivatives accounting for quaternion specialness
# @inline function ∂g∂ʳvela(joint::Rotational{T,N}, q1a::UnitQuaternion, xa::AbstractVector,
#     qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
#     q1b::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion,
#     vb::AbstractVector, ωb::AbstractVector, Δt) where {T,N}
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#     τ_spring = springtorque(joint, qa, qb)
#     τ_damp = dampertorque(joint, qa, ωa, qb, ωb)
#     qoffset = joint.qoffset
#     Vdamp = szeros(T, 3, 3)
#     Qdamp = ∂vrotate∂p(τ_damp, qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(qb * inv(qoffset)) * Tmat()
#     Qdamp += ∂vrotate∂q(τ_damp, qa * qoffset) * Rmat(qoffset)
#     Ωdamp = Qdamp * Lmat(q1a) * derivωbar(ωa, Δt) * Δt/2
#     Ωdamp += ∂vrotate∂p(τ_damp, qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * -1.0 * ∂vrotate∂p(ωa, inv(qoffset))
#
#     Vspring = szeros(T, 3, 3)
#     Qspring = ∂vrotate∂p(τ_spring, qa * qoffset) * -Aᵀ * A * joint.spring * Aᵀ * A * VRmat(qb * inv(qoffset)) * Tmat()
#     Qspring += ∂vrotate∂q(τ_spring, qa * qoffset) * Rmat(qoffset)
#     Ωspring = Qspring * Lmat(q1a) * derivωbar(ωa, Δt) * Δt/2
#     # Ωspring += nothing
#
#     V = Vspring + Vdamp
#     Ω = Ωspring + Ωdamp
#     V *= Δt
#     Ω *= Δt
#
#     return [Aᵀ * A * V Aᵀ * A * Ω]
# end
# @inline function ∂g∂ʳvelb(joint::Rotational{T,N}, q1a::UnitQuaternion, xa::AbstractVector,
#     qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
#     q1b::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion,
#     vb::AbstractVector, ωb::AbstractVector, Δt) where {T,N}
#
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#     τ_spring = springtorque(joint, qa, qb)
#     τ_damp = dampertorque(joint, qa, ωa, qb, ωb)
#     qoffset = joint.qoffset
#     Vdamp = szeros(T, 3, 3)
#     Qdamp = ∂vrotate∂p(τ_damp, qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(inv(qoffset)) * Lmat(inv(qa))
#     Ωdamp = Qdamp * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2
#     Ωdamp += ∂vrotate∂p(τ_damp, qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂p(ωb, qa \ qb / qoffset)
#
#     Vspring = szeros(T, 3, 3)
#     Qspring = ∂vrotate∂p(τ_spring, qa * qoffset) * -Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * Lmat(inv(qa))
#     Ωspring = Qspring * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2
#     # Ωspring += nothing
#
#     V = Vspring + Vdamp
#     Ω = Ωspring + Ωdamp
#     V *= Δt
#     Ω *= Δt
#
#     return [Aᵀ * A * V Aᵀ * A * Ω]
# end
# @inline function ∂g∂ʳvelb(joint::Rotational{T,N}, q1b::UnitQuaternion, xb::AbstractVector,
# qb::UnitQuaternion, vb::AbstractVector,  ωb::AbstractVector, Δt) where {T,N}
#         A = nullspacemat(joint)
#         Aᵀ = zerodimstaticadjoint(A)
#         τ_spring = springtorque(joint, qb)
#         τ_damp = dampertorque(joint, qb, ωb)
#         qoffset = joint.qoffset
#         Vdamp = szeros(T, 3, 3)
#         Qdamp = ∂vrotate∂p(τ_damp, qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qb / qoffset) * Rmat(inv(qoffset))
#         Ωdamp = Qdamp * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2
#         Ωdamp += ∂vrotate∂p(τ_damp, qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂p(ωb, qb / qoffset)
#
#         Vspring = szeros(T, 3, 3)
#         Qspring = -1.0 * ∂vrotate∂p(τ_spring, qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset))
#         Ωspring = Qspring * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2
#         # Ωspring += nothing
#
#         V = Vspring + Vdamp
#         Ω = Ωspring + Ωdamp
#         V *= Δt
#         Ω *= Δt
#     return [Aᵀ * A * V Aᵀ * A * Ω]
# end
#
# ## vec(G) Jacobian (also NOT accounting for quaternion specialness in the second derivative: ∂(∂ʳg∂posx)∂y)
# @inline function ∂2g∂posaa(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
#     Lpos = Lmat(UnitQuaternion(xb + vrotate(joint.vertices[2], qb) - xa))
#     Ltpos = Lᵀmat(UnitQuaternion(xb + vrotate(joint.vertices[2], qb) - xa))
#
#     XX = szeros(T, 9, 3) # empty
#     XQ = szeros(T, 9, 4) # empty
#     QX = szeros(T, 9, 3) # empty
#
#     f = q -> ∂g∂ʳposa(joint, xa, UnitQuaternion(q...), xb, qb)[1:3, 4:6]
#     df = ForwardDiff.jacobian(f, [qa.w; qa.x; qa.y; qa.z])
#     # @show df
#
#     QQ = df#szeros(T, 9, 4)
#
#     return XX, XQ, QX, QQ
# end
# @inline function ∂2g∂posab(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
#     XX = szeros(T, 9, 3) # empty
#     XQ = szeros(T, 9, 4) # empty
#     QX = szeros(T, 9, 3) # empty
#
#     f = q -> ∂g∂ʳposa(joint, xa, qa, xb, UnitQuaternion(q...))[1:3, 4:6]
#     df = ForwardDiff.jacobian(f, [qb.w; qb.x; qb.y; qb.z])
#     # @show df
#     QQ = df#szeros(T, 9, 4)
#
#     return XX, XQ, QX, QQ
# end
# @inline function ∂2g∂posba(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
#     XX = szeros(T, 9, 3) # empty
#     XQ = szeros(T, 9, 4) # empty
#     QX = szeros(T, 9, 3) # empty
#
#     f = q -> ∂g∂ʳposb(joint, xa, UnitQuaternion(q...), xb, qb)[1:3, 4:6]
#     df = ForwardDiff.jacobian(f, [qa.w; qa.x; qa.y; qa.z])
#     # @show df
#
#     QQ = df#szeros(T, 9, 4)
#
#     return XX, XQ, QX, QQ
# end
# @inline function ∂2g∂posbb(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
#     XX = szeros(T, 9, 3) # empty
#     XQ = szeros(T, 9, 4) # empty
#     QX = szeros(T, 9, 3) #empty
#
#     f = q -> ∂g∂ʳposb(joint, xa, qa, xb, UnitQuaternion(q...))[1:3, 4:6]
#     df = ForwardDiff.jacobian(f, [qb.w; qb.x; qb.y; qb.z])
#     # @show df
#
#     QQ = df#szeros(T, 9, 4)
#
#
#     return XX, XQ, QX, QQ
# end
# @inline function ∂2g∂posbb(joint::Rotational{T}, xb::AbstractVector, qb::UnitQuaternion) where T
#     XX = szeros(T, 9, 3) # empty
#     XQ = szeros(T, 9, 4) # empty
#     QX = szeros(T, 9, 3) # empty
#
#     f = q -> ∂g∂ʳposb(joint, xb, UnitQuaternion(q...))[1:3, 4:6]
#     df = ForwardDiff.jacobian(f, [qb.w; qb.x; qb.y; qb.z])
#
#     # @show df
#     QQ = df#szeros(T, 9, 4)
#
#     return XX, XQ, QX, QQ
# end
#
# ### Forcing
# ## Application of joint forces (for dynamics)
# @inline function applyFτ!(joint::Rotational{T}, statea::State, stateb::State, Δt::T, clear::Bool) where T
#     τ = joint.Fτ
#     _, qa = posargsk(statea)
#     _, qb = posargsk(stateb)
#
#     τa = -τ # in world coordinates
#     τb = τ  # in world coordinates
#
#     τa = vrotate(τa, inv(qa)) # in local coordinates
#     τb = vrotate(τb, inv(qb)) # in local coordinates
#
#     statea.τk[end] += τa
#     stateb.τk[end] += τb
#     clear && (joint.Fτ = szeros(T,3))
#     return
# end
# @inline function applyFτ!(joint::Rotational{T}, stateb::State, Δt::T, clear::Bool) where T
#     τ = joint.Fτ
#     _, qb = posargsk(stateb)
#
#     τb = τ # in world coordinates
#
#     τb = vrotate(τb,inv(qb)) # in local coordinates
#
#     stateb.τk[end] += τb
#     clear && (joint.Fτ = szeros(T,3))
#     return
# end
#
# ## Forcing derivatives (for linearization)
# # Control derivatives
# @inline function ∂Fτ∂ua(joint::Rotational, statea::State, stateb::State)
#     BFa = szeros(T, 3, 3)
#     Bτa = szeros(T, 3, 3)
#
#     return [BFa; Bτa]
# end
#
# @inline function ∂Fτ∂ub(joint::Rotational, statea::State, stateb::State)
#     BFb = szeros(T, 3, 3)
#     Bτb = szeros(T, 3, 3)
#
#     return [BFb; Bτb]
# end
# @inline function ∂Fτ∂ub(joint::Rotational, stateb::State)
#     BFb = szeros(T, 3, 3)
#     Bτb = szeros(T, 3, 3)
#     return [BFb; Bτb]
# end
#
# # Position derivatives
# @inline function ∂Fτ∂posa(joint::Rotational{T}, statea::State, stateb::State) where T
#     FaXa = szeros(T,3,3)
#     FaQa = szeros(T,3,3)
#     τaXa = szeros(T,3,3)
#     τaQa = szeros(T,3,3)
#     FbXa = szeros(T,3,3)
#     FbQa = szeros(T,3,3)
#     τbXa = szeros(T,3,3)
#     τbQa = szeros(T,3,3)
#
#     return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
# end
# @inline function ∂Fτ∂posb(joint::Rotational{T}, statea::State, stateb::State) where T
#     FaXb = szeros(T,3,3)
#     FaQb = szeros(T,3,3)
#     τaXb = szeros(T,3,3)
#     τaQb = szeros(T,3,3)
#     FbXb = szeros(T,3,3)
#     FbQb = szeros(T,3,3)
#     τbXb = szeros(T,3,3)
#     τbQb = szeros(T,3,3)
#
#     return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
# end
# @inline function ∂Fτ∂posb(joint::Rotational{T}, stateb::State) where T
#     xb, qb = posargsk(stateb)
#     F = joint.Fτ
#     vertices = joint.vertices
#
#     FaXb = szeros(T,3,3)
#     FaQb = szeros(T,3,3)
#     τaXb = szeros(T,3,3)
#     τaQb = szeros(T,3,3)
#     FbXb = szeros(T,3,3)
#     FbQb = szeros(T,3,3)
#     τbXb = szeros(T,3,3)
#     τbQb = szeros(T,3,3)
#
#     return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
# end
#
#
# # ### Minimal coordinates
# # ## Position and velocity offsets
# # @inline function getPositionDelta(joint::Rotational, body1::AbstractBody, body2::Body, x::SVector)
# #     Δx = zerodimstaticadjoint(nullspacemat(joint)) * x # in body1 frame
# #     return Δx
# # end
# # @inline function getVelocityDelta(joint::Rotational, body1::AbstractBody, body2::Body, v::SVector)
# #     Δv = zerodimstaticadjoint(nullspacemat(joint)) * v # in body1 frame
# #     return Δv
# # end
#
# # ## Minimal coordinate calculation
# # @inline function minimalCoordinates(joint::Rotational, body1::Body, body2::Body)
# #     statea = body1.state
# #     stateb = body2.state
# #     return nullspacemat(joint) * g(joint, statea.xc, statea.qc, stateb.xc, stateb.qc)
# # end
# # @inline function minimalCoordinates(joint::Rotational, body1::Origin, body2::Body)
# #     stateb = body2.state
# #     return nullspacemat(joint) * g(joint, stateb.xc, stateb.qc)
# # end
# # @inline function minimalVelocities(joint::Rotational, body1::Body, body2::Body)
# #     statea = body1.state
# #     stateb = body2.state
# #     return nullspacemat(joint) * (stateb.vc - statea.vc)
# # end
# # @inline function minimalVelocities(joint::Rotational, body1::Origin, body2::Body)
# #     stateb = body2.state
# #     return nullspacemat(joint) * stateb.vc
# # end
# # nullspacemat(force1)
#
# function ∂springtorquea∂posa(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
#
# end
#
# function ∂dampertorquea∂posa(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
#
# end
#
# function ∂springtorqueb∂posb(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
#
# end
#
# function ∂dampertorqueb∂posb(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
#
# end
#
# function ∂springtorqueb∂posb(joint::Rotational, body1::Origin, body2::Body, Δt::T) where T
#
# end
#
# function ∂dampertorqueb∂posb(joint::Rotational, body1::Origin, body2::Body, Δt::T) where T
#
# end
