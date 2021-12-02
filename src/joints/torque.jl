#
## Position level constraints (for dynamics) in world frame
"""
    torque = 1/2 * spring * sin(θ) ≡≡ k * w * sqrt(1 - w^2)
    q = [sin(θ/2); r cos(θ/2)] = [w, x, y, z]
"""
@inline function gc(joint::Rotational, qa::UnitQuaternion, qb::UnitQuaternion; qoff::UnitQuaternion = one(UnitQuaternion))
    # return Vmat(qa \ qb / joint.qoffset)

    # q = qa \ qb / joint.qoffset
    # angle = 2 * acos(q.w)
    # axis = Vmat(q) ./ sqrt(1 - q.w*q.w)
    # return angle * axis

    q = qa \ qb / joint.qoffset# / qoff
    return Vmat(q) * q.w
end
"""
    torque = 1/2 * spring * sin(θ) ≡≡ k * w * sqrt(1 - w^2)
    q = [sin(θ/2); r cos(θ/2)] = [w, x, y, z]
"""
@inline function gc(joint::Rotational, qb::UnitQuaternion; qoff::UnitQuaternion = one(UnitQuaternion))
    # return Vmat(qb / joint.qoffset)

    # q = qb / joint.qoffset
    # angle = 2 * acos(q.w)
    # axis = Vmat(q) ./ sqrt(1 - q.w*q.w)
    # return angle * axis

    q = qb / joint.qoffset# / qoff
    return Vmat(q) * q.w
end
# used to compute potential energy
@inline function gq(joint::Rotational, qa::UnitQuaternion, qb::UnitQuaternion)
    return qa \ qb / joint.qoffset
end
@inline function gq(joint::Rotational, qb::UnitQuaternion)
    return qb / joint.qoffset
end

### Spring and damper
## Discrete-time position wrappers (for dynamics)
springforcea(joint::Rotational, statea::State, stateb::State, Δt) = Δt * springforcea(joint, posargs2(statea)[2], posargs2(stateb)[2])
springforceb(joint::Rotational, statea::State, stateb::State, Δt) = Δt * springforceb(joint, posargs2(statea)[2], posargs2(stateb)[2])
springforceb(joint::Rotational, stateb::State, Δt) = Δt * springforceb(joint, posargs2(stateb)[2])
damperforcea(joint::Rotational, statea::State, stateb::State, Δt) = Δt * damperforcea(joint, posargs2(statea)[2], statea.ϕsol[2], posargs2(stateb)[2], stateb.ϕsol[2])
damperforceb(joint::Rotational, statea::State, stateb::State, Δt) = Δt * damperforceb(joint, posargs2(statea)[2], statea.ϕsol[2], posargs2(stateb)[2], stateb.ϕsol[2])
damperforceb(joint::Rotational, stateb::State, Δt) = Δt * damperforceb(joint, posargs2(stateb)[2], stateb.ϕsol[2])

springforcea(joint::Rotational{T,3}, statea::State, stateb::State, Δt) where {T} = szeros(T, 6)
springforceb(joint::Rotational{T,3}, statea::State, stateb::State, Δt) where {T} = szeros(T, 6)
springforceb(joint::Rotational{T,3}, stateb::State, Δt) where {T} = szeros(T, 6)
damperforcea(joint::Rotational{T,3}, statea::State, stateb::State, Δt) where {T} = szeros(T, 6)
damperforceb(joint::Rotational{T,3}, statea::State, stateb::State, Δt) where {T} = szeros(T, 6)
damperforceb(joint::Rotational{T,3}, stateb::State, Δt) where {T} = szeros(T, 6)

# Used in energy computation
springforcea(joint::Rotational{T,3}, qa::UnitQuaternion, qb::UnitQuaternion; rotate::Bool = true) where {T} = szeros(T, 6)
springforceb(joint::Rotational{T,3}, qa::UnitQuaternion, qb::UnitQuaternion; rotate::Bool = true) where {T} = szeros(T, 6)
springforceb(joint::Rotational{T,3}, qb::UnitQuaternion; rotate::Bool = true) where {T} = szeros(T, 6)

### Spring and damper
# Force applied by body b on body a expressed in frame a
@inline function springforcea(joint::Rotational{T}, qa::UnitQuaternion, qb::UnitQuaternion; rotate::Bool = true) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    # We need to convert the joint.spring_offset which is some minimal coordinate representation of the joint. For us it's the axis angle reresentation.
    # We convert it to the 'sinusoidal spring'
    aa = Aᵀ * joint.spring_offset # axis angle
    θ = norm(aa)
    qoff = UnitQuaternion(cos(θ/2), 1/2 * sinc(θ/(2π)) * aa) # quaternion
    # offset = Vmat(qoff) * qoff.w
    # distance = A * (gc(joint, qa, qb) .- offset)
    # distance = A * gc(joint, qa, qb, qoff = qoff)
    distance = A * gc(joint, qa, qb)
    force = joint.spring * Aᵀ * distance # force in offset frame
    rotate && (force = vrotate(force, joint.qoffset)) # rotate back to a frame
    return [szeros(T, 3); force]
end
# Force applied by body a on body b expressed in frame b
@inline function springforceb(joint::Rotational{T}, qa::UnitQuaternion, qb::UnitQuaternion; rotate::Bool = true) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    # We need to convert the joint.spring_offset which is some minimal coordinate representation of the joint. For us it's the axis angle reresentation.
    # We convert it to the 'sinusoidal spring'
    aa = Aᵀ * joint.spring_offset # axis angle
    θ = norm(aa)
    qoff = UnitQuaternion(cos(θ/2), 1/2 * sinc(θ/(2π)) * aa) # quaternion
    # offset = Vmat(qoff) * qoff.w
    # distance = A * (gc(joint, qa, qb) .- offset)
    # distance = A * gc(joint, qa, qb, qoff = qoff)
    distance = A * gc(joint, qa, qb)
    force = - joint.spring * Aᵀ * distance # force in offset frame
    rotate && (force = vrotate(force, inv(qb) * qa * joint.qoffset)) # rotate back to b frame
    return [szeros(T, 3); force]
end
# Force applied by origin on body b expressed in frame b
@inline function springforceb(joint::Rotational{T}, qb::UnitQuaternion; rotate::Bool = true) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    # We need to convert the joint.spring_offset which is some minimal coordinate representation of the joint. For us it's the axis angle reresentation.
    # We convert it to the 'sinusoidal spring'
    aa = Aᵀ * joint.spring_offset # axis angle
    θ = norm(aa)
    qoff = UnitQuaternion(cos(θ/2), 1/2 * sinc(θ/(2π)) * aa) # quaternion
    # offset = Vmat(qoff) * qoff.w
    # distance = A * (gc(joint, qb) .- offset)
    # distance = A * gc(joint, qb, qoff = qoff)
    distance = A * gc(joint, qb)
    force = - joint.spring * Aᵀ * distance # force in offset frame
    rotate && (force = vrotate(force, inv(qb) * joint.qoffset)) # rotate back to b frame
    return [szeros(T, 3); force]
end

# Force applied by body b on body a expressed in frame a
@inline function damperforcea(joint::Rotational{T}, qa::UnitQuaternion, ωa::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector; rotate::Bool = true) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    qoffset = joint.qoffset
    velocity = A * (vrotate(ωb, qa \ qb / qoffset) - vrotate(ωa, inv(qoffset))) # in offset frame
    force = 2 * Aᵀ * A * joint.damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
    rotate && (force = vrotate(force, qoffset)) # rotate back to frame a
    return [szeros(T, 3); force]
end
# Force applied by body a on body b expressed in frame b
@inline function damperforceb(joint::Rotational{T}, qa::UnitQuaternion, ωa::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector; rotate::Bool = true) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    qoffset = joint.qoffset
    velocity = A * (vrotate(ωb, qa \ qb / qoffset) - vrotate(ωa, inv(qoffset))) # in offset frame
    force = - 2 * Aᵀ * A * joint.damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
    rotate && (force = vrotate(force, inv(qb) * qa * qoffset)) # rotate back to frame b
    return [szeros(T, 3); force]
end
# Force applied by origin on body b expressed in frame b
@inline function damperforceb(joint::Rotational{T}, qb::UnitQuaternion, ωb::AbstractVector; rotate::Bool = true) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    qoffset = joint.qoffset
    velocity = A * vrotate(ωb, qb / qoffset) # in offset frame
    force = - 2 * Aᵀ * A * joint.damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
    rotate && (force = vrotate(force, inv(qb) * qoffset)) # rotate back to frame b
    return [szeros(T, 3); force]
end

∂springforcea∂posa(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂damperforcea∂posa(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂springforcea∂posb(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂damperforcea∂posb(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂springforceb∂posb(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂damperforceb∂posb(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂springforceb∂posa(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂damperforceb∂posa(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂springforceb∂posb(joint::Rotational{T,3}, body1::Origin, body2::Body, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂damperforceb∂posb(joint::Rotational{T,3}, body1::Origin, body2::Body, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)

∂springforcea∂vela(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T) where T = szeros(T, 6, 6)
∂damperforcea∂vela(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T) where T = szeros(T, 6, 6)
∂springforcea∂velb(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T) where T = szeros(T, 6, 6)
∂damperforcea∂velb(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T) where T = szeros(T, 6, 6)
∂springforceb∂velb(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T) where T = szeros(T, 6, 6)
∂damperforceb∂velb(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T) where T = szeros(T, 6, 6)
∂springforceb∂vela(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T) where T = szeros(T, 6, 6)
∂damperforceb∂vela(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T) where T = szeros(T, 6, 6)
∂springforceb∂velb(joint::Rotational{T,3}, body1::Origin, body2::Body, Δt::T) where T = szeros(T, 6, 6)
∂damperforceb∂velb(joint::Rotational{T,3}, body1::Origin, body2::Body, Δt::T) where T = szeros(T, 6, 6)

function ∂springforcea∂posa(joint::Rotational, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    # qoffset = joint.qoffset
    # force = springforcea(joint, qa, qb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    # Q = ∂vrotate∂p(force, qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(qb * inv(qoffset)) * Tmat() * LVᵀmat(qa)
    Q = FiniteDiff.finite_difference_jacobian(qa -> springforcea(joint, UnitQuaternion(qa..., false), qb)[SVector{3,Int}(4,5,6)], [qa.w, qa.x, qa.y, qa.z])
    attjac && (Q *= LVᵀmat(qa))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end
function ∂damperforcea∂posa(joint::Rotational, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforcea(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, qoffset) * 2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(qb * inv(qoffset)) * Tmat()
    attjac && (Q *= LVᵀmat(qa))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end
function ∂springforcea∂posb(joint::Rotational, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    # qoffset = joint.qoffset
    # force = springforcea(joint, qa, qb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    # Q = ∂vrotate∂p(force, qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * Lmat(inv(qa)) * LVᵀmat(qb)
    Q = FiniteDiff.finite_difference_jacobian(qb -> springforcea(joint, qa, UnitQuaternion(qb..., false))[SVector{3,Int}(4,5,6)], [qb.w, qb.x, qb.y, qb.z])
    attjac && (Q *= LVᵀmat(qb))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end
function ∂damperforcea∂posb(joint::Rotational, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1b, q1b = posargs2(body2.state)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforcea(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, qoffset) * 2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(inv(qoffset)) * Lmat(inv(qa))
    attjac && (Q *= LVᵀmat(qb))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end
function ∂springforceb∂posb(joint::Rotational, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    # qoffset = joint.qoffset
    # force = springforceb(joint, qa, qb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    # Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -1.0 * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * Lmat(inv(qa)) * LVᵀmat(qb)
    # Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qa * qoffset) * Tmat() * LVᵀmat(qb)
    Q = FiniteDiff.finite_difference_jacobian(qb -> springforceb(joint, qa, UnitQuaternion(qb..., false))[SVector{3,Int}(4,5,6)], [qb.w, qb.x, qb.y, qb.z])
    attjac && (Q *= LVᵀmat(qb))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end
function ∂damperforceb∂posb(joint::Rotational, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforceb(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(inv(qoffset)) * Lmat(inv(qa))
    Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qa * qoffset) * Tmat()
    attjac && (Q *= LVᵀmat(qb))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end
function ∂springforceb∂posa(joint::Rotational, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # _, _, _, ωa = fullargssol(body1.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    # qoffset = joint.qoffset
    # force = springforceb(joint, qa, qb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    # Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -1.0 * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(qb * inv(qoffset)) * Tmat() * LVᵀmat(qa)
    # Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qoffset) * Lmat(inv(qb)) * LVᵀmat(qa)
    Q = FiniteDiff.finite_difference_jacobian(qa -> springforceb(joint, UnitQuaternion(qa..., false), qb)[SVector{3,Int}(4,5,6)], [qa.w, qa.x, qa.y, qa.z])
    attjac && (Q *= LVᵀmat(qa))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end
function ∂damperforceb∂posa(joint::Rotational, body1::Body, body2::Body, Δt::T; attjac::Bool = true) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforceb(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(qb * inv(qoffset)) * Tmat()
    Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qoffset) * Lmat(inv(qb))
    attjac && (Q *= LVᵀmat(qa))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end
function ∂springforceb∂posb(joint::Rotational, body1::Origin, body2::Body, Δt::T; attjac::Bool = true) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # _, _, _, ωb = fullargssol(body2.state)
    xb, qb = posargs2(body2.state)
    # qoffset = joint.qoffset
    # force = springforceb(joint, qb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    # Q = ∂vrotate∂p(force, inv(qb) * qoffset) * -1.0 * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * LVᵀmat(qb)
    # Q += ∂vrotate∂q(force, inv(qb) * qoffset) * Rmat(qoffset) * Tmat() * LVᵀmat(qb)
    Q = FiniteDiff.finite_difference_jacobian(qb -> springforceb(joint, UnitQuaternion(qb..., false))[SVector{3,Int}(4,5,6)], [qb.w, qb.x, qb.y, qb.z])
    attjac && (Q *= LVᵀmat(qb))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end
function ∂damperforceb∂posb(joint::Rotational, body1::Origin, body2::Body, Δt::T; attjac::Bool = true) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωb = fullargssol(body2.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforceb(joint, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, inv(qb) * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qb / qoffset) * Rmat(inv(qoffset))
    Q += ∂vrotate∂q(force, inv(qb) * qoffset) * Rmat(qoffset) * Tmat()
    attjac && (Q *= LVᵀmat(qb))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end


function ∂springforcea∂vela(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    V = szeros(T, 3, 3)
    Ω = szeros(T, 3, 3)
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂damperforcea∂vela(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforcea(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, qoffset) * 2 * Aᵀ * A * joint.damper * Aᵀ * A * -1.0 * ∂vrotate∂p(ωa, inv(qoffset))
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂springforcea∂velb(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    V = szeros(T, 3, 3)
    Ω = szeros(T, 3, 3)
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂damperforcea∂velb(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforcea(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, qoffset) * 2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂p(ωb, qa \ qb / qoffset)
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂springforceb∂velb(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    V = szeros(T, 3, 3)
    Ω = szeros(T, 3, 3)
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂damperforceb∂velb(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforceb(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂p(ωb, qa \ qb / qoffset)
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂springforceb∂vela(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    V = szeros(T, 3, 3)
    Ω = szeros(T, 3, 3)
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂damperforceb∂vela(joint::Rotational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforceb(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * -1.0 * ∂vrotate∂p(ωa, inv(qoffset))
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂springforceb∂velb(joint::Rotational, body1::Origin, body2::Body, Δt::T) where T
    V = szeros(T, 3, 3)
    Ω = szeros(T, 3, 3)
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂damperforceb∂velb(joint::Rotational, body1::Origin, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωb = fullargssol(body2.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforceb(joint, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, inv(qb) * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂p(ωb, qb / qoffset)
    return Δt * [szeros(T, 3, 6); V Ω]
end
